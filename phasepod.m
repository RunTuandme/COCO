function [final_orbit, final_clk, final_amb, orbit_list, N_LOOP, X_correct, RMSE, inter_t, y_error] = phasepod(phasedata, pseudodata, init_orbit, JD2_init, t0_init, ts)
% DESCRIPTION:     Presion Orbit Determination.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT：          phasedata pseudodata (format)
%                 | phase  | epoch | Ambiguity Number | GNSS number | GNSS x | GNSS y | GNSS z | GNSS clk | receiver clk | NaN | NaN | Second(INT) | JD1 | JD2 |
%                 | pseudo | epoch | NaN              | GNSS number | GNSS x | GNSS y | GNSS z | GNSS clk | receiver clk | NaN | NaN | Second(INT) | JD1 | JD2 |
%                 init_orbit, Initial orbit
%                 t_init, Julian date of initial orbit
%                 t0_init, Second(INT) of initial orbit
%                 ts, interval of data

%% 初始化并载入全局变量
c = 2.99792458E+8;           % m/s
omega_earth = 7.292115146706979e-5; % rad/s

N_phasedata = size(phasedata, 1);
N_pseudodata = size(pseudodata, 1);
N_all = N_phasedata + N_pseudodata;

phasedata(:,10) = (phasedata(:,2) - 1) * ts - t0_init;
phasedata(:,11) = 99 * N_all / (99 * N_phasedata + N_pseudodata);
pseudodata(:,10) = (pseudodata(:,2) - 1) * ts - t0_init;
pseudodata(:,11) = N_all / (99 * N_phasedata + N_pseudodata);

RMSE = Inf;
RMSE_new = 0;

amb_num = unique(phasedata(:,3));    % 模糊度个数
N_phase_ambiguity = numel(amb_num);

combinedata = [phasedata; pseudodata];
combinedata(:,15) = combinedata(:,3);
% 将第15列中的所有NaN值替换为-1
column_15 = combinedata(:, 15);           
nan_indices = isnan(column_15);       
column_15(nan_indices) = -1;            
combinedata(:, 15) = column_15;            
combinedata = sortrows(combinedata, [2,4,15]); %按历元排序

combinedata(:,16) = 1;

% combinedata = combinedata(mod(combinedata(:, 2), 10) == 1, :);
% combinedata = combinedata(combinedata(:, 2) < 21600, :);

t_end = combinedata(end,10);
clk_num = unique(combinedata(:,2));    % 钟差个数，对应每个历元
N_clk_error = numel(clk_num);

% CdAM = 0.044;
% CrAM = 0.02;
CdAM = 0.01;
CrAM = 0.012;

% 大气光压参数分段
global Cd_tinter
global Cr_tinter
global num_Cd
global num_Cr
global exp_switch
Cd_tinter = 3600 * 3;  % second，不估可以设 inf
Cr_tinter = 3600 * 12; % second
num_Cd = ceil(t_end / Cd_tinter);
num_Cr = ceil(t_end / Cr_tinter);

CC = zeros(3,1);
Cn = zeros(3,1);
Sn = zeros(3,1);

if exp_switch
    N_para = num_Cd + num_Cr + 9; % 大气1、光压1、经验力9
else
    N_para = num_Cd + num_Cr;
    CC = [];
    Cn = [];
    Sn = [];
end
X_correct = zeros(6+N_para,1);  % 状态量
H_init = [eye(6) zeros(6,N_para)]; 
H_init = [X_correct; H_init(:)]; 

amb_state = zeros(N_phase_ambiguity, 1);
clk_state = zeros(N_clk_error, 1);
%钟差初始化
for i = 1:N_clk_error
    clki = combinedata(combinedata(:,2) == clk_num(i),9);
    clk_state(i) = clki(1);
end

init_state = [init_orbit; CdAM*ones(num_Cd,1); CrAM*ones(num_Cr,1); CC; Cn; Sn; amb_state; clk_state];   % 当前估计值
init_orbit = [init_orbit; CdAM*ones(num_Cd,1); CrAM*ones(num_Cr,1); CC; Cn; Sn; amb_state; clk_state];   % 先验值

opts = odeset('Reltol',1e-12,'AbsTol',1e-13,'Stats','off','NormControl','on');

% TAU = 0.001239447665903; % 1/s, sqrt(GM/ae^3)
% LEN = 6.37813630E+06;    % m,   ae
% scale_STM = [ones(3,3), TAU * ones(3,3), ones(3,N_para) / LEN;
%              ones(3,3) / TAU, ones(3,3), ones(3,N_para) / LEN / TAU;
%              ones(N_para,3) * LEN, ones(N_para,3) * LEN * TAU, ones(N_para,N_para)];


N_LOOP = 0;
ode_flag = true;
residual_limit = 0;
sigma3_p = Inf;
sigma3_l = Inf;
%% 计算轨道积分&状态转移矩阵
while((residual_limit < 2) && (N_LOOP < 20))
    N_LOOP = N_LOOP + 1;
    fprintf('开始第%d次迭代...\n', N_LOOP);
    
    % 更新phasedata和pseudodata
    N_all = size(combinedata, 1);
    if N_phasedata + N_pseudodata ~= N_all
        phasedata = combinedata(~isnan(combinedata(:,3)),:);
        pseudodata = combinedata(isnan(combinedata(:,3)),:);
    end
    
    % 剔除后删减钟差、模糊度参数
    N_phasedata = size(phasedata, 1);
    N_pseudodata = size(pseudodata, 1);
    if N_clk_error ~= numel(unique(combinedata(:,2)))
        amb_state_before = init_state(6+N_para+1:6+N_para+N_phase_ambiguity);
        clk_state_before = init_state(6+N_para+N_phase_ambiguity+1:end);
        
        clk_num_new = unique(combinedata(:,2));
        [~, clk_pos] = ismember(clk_num_new, clk_num); % 获取clk_num_new中元素在clk_num中的位置
        
        clk_state = clk_state_before(clk_pos);
        init_state = [init_state(1:6+N_para); amb_state_before; clk_state];
        
        N_clk_error = numel(unique(combinedata(:,2)));
        clk_num = clk_num_new;
    end
    if N_phase_ambiguity ~= numel(unique(phasedata(:,3)))
        amb_state_before = init_state(6+N_para+1:6+N_para+N_phase_ambiguity);
        clk_state_before = init_state(6+N_para+N_phase_ambiguity+1:end);
        
        amb_num_new = unique(phasedata(:,3));
        [~, amb_pos] = ismember(amb_num_new, amb_num);
        
        amb_state = amb_state_before(amb_pos);
        init_state = [init_state(1:6+N_para); amb_state; clk_state_before];
        
        N_phase_ambiguity = numel(unique(phasedata(:,3)));
        amb_num = amb_num_new;
    end
    
    phasedata(:,11) = 99 * N_all / (99 * N_phasedata + N_pseudodata);
    pseudodata(:,11) = N_all / (99 * N_phasedata + N_pseudodata);
    combinedata(~isnan(combinedata(:,3)),11) = 99 * N_all / (99 * N_phasedata + N_pseudodata);
    combinedata(isnan(combinedata(:,3)),11) = N_all / (99 * N_phasedata + N_pseudodata);

    if ode_flag
        tic
        disp('轨道积分...');
        % 全局参数
        H_init(1:6+N_para) = init_state(1:6+N_para);

        [time, vector_y] = ode113(@variationalEquation, [0 t_end], H_init, opts, JD2_init);
        toc
    end
    
    inter_t = unique(combinedata(:,10)); % 第 10 列为积分时长，根据历元号对应得到
%     vector_yin = interp1(time, vector_y, inter_t, 'spline');

    N_inter = length(inter_t);
    N_ycol = size(vector_y, 2);
    vector_yin = zeros(N_inter, N_ycol);
    for i = 1:N_ycol
        ilist = interp_larange_list(time, vector_y(:,i), inter_t);
        vector_yin(:,i) = ilist;
    end

   %% 计算加权法方程
    y_c = zeros(N_all, 1);
    ri = zeros(N_all, 1);
    %HGi = zeros(3, 3, N_all);
    rrhoi = zeros(N_all ,3);
    disp("计算残差")
    tic
    % 计算残差
    clk_list = unique(combinedata(:,2));
    amb_list = unique(phasedata(:,3));
    amb_state = init_state(6+N_para+1 : 6+N_para+N_phase_ambiguity);
    clk_state = init_state(end-N_clk_error+1 : end);

    for i = 1: N_all
        % get value from combinedata
        rho_obs = combinedata(i,1);
        epoch = combinedata(i,2);
        amb = combinedata(i,3);
        gnssnum = combinedata(i,4);
        ti = combinedata(i,12);
        JD1 = combinedata(i,13);
        JD2 = combinedata(i,14);
  
        dataN = find(clk_list==epoch);
        r_leopre_GCRS = vector_yin(dataN, 1:3)';
        v_leopre_GCRS = vector_yin(dataN, 4:6)';
        clk_leo = clk_state(dataN);
        if isnan(amb)
            ambdata = 0;
        else
            ambdata = amb_state(amb_list==amb);
        end
        
        HG = sofa_C2T(JD1, JD2);
        r_leopre_ITRS = HG * r_leopre_GCRS;
        [f3d, satfixed2earthfixed] = TOOL_getapcoffset(JD1, JD2, r_leopre_ITRS, HG);
        r_leopre_ITRS = r_leopre_ITRS + f3d;
        r_leopre_GCRS = HG' * r_leopre_ITRS;
             
        if i ~= 1
            if gnssnum == combinedata(i-1,4) && epoch == combinedata(i-1,2)
                % Do nothing
                % rd and r_gnss_ITRS_row retain their values from the last loop
            else
                [rd, r_gnss_ITRS_row, r_rho] = TOOL_lighttime(ti, JD1, JD2, gnssnum, r_leopre_GCRS, v_leopre_GCRS, clk_leo, HG);
            end
        else
            [rd, r_gnss_ITRS_row, r_rho] = TOOL_lighttime(ti, JD1, JD2, gnssnum, r_leopre_GCRS, v_leopre_GCRS, clk_leo, HG);
        end
        
        POS_satfixed = satfixed2earthfixed' * (r_gnss_ITRS_row'-r_leopre_ITRS);
        dxy = sqrt((POS_satfixed(2))^2 + (POS_satfixed(1))^2);
        angel = atan2(POS_satfixed(3), dxy) * 180 / pi;
        
        combinedata(i,5:7) = r_gnss_ITRS_row;
        combinedata(i,17) = angel;
        %HGi(:, :, i) = HG;
        rrhoi(i, :) = r_rho'; 
        ri(i) = rd;
        y_c(i) = rho_obs - (rd - c * TOOL_getCLK(ti, gnssnum) + clk_leo + ambdata);
    end
    toc
    
    % 计算法方程
    tic
    disp('计算法方程');
    fepo = 0;
    clk_mi = 0;
    B = 0;
    L = 0;
    H0yi = cell(N_clk_error, 3);
    for i = 1:N_all
        epoch = combinedata(i,2);
        % 历元开始
        if epoch ~= fepo || i == N_all
            if i == N_all
                dYdamb = zeros(1, N_phase_ambiguity);
                if ~isnan(combinedata(i,3))
                    dYdamb(amb_list == combinedata(i,3)) = 1;
                end
                dYdXnew = [rrhoi(i,:) / ri(i), ...
                    0, 0, 0, ...
                    zeros(1, N_para), ...
                    dYdamb];
                dYdX = [dYdX; dYdXnew];
                P = [P; combinedata(i, 11)];
                i = i + 1;
            end
            % 计算上一历元B、L矩阵
            if fepo
                y = vector_yin(clk_mi,:)';
                STM = [reshape(y(6+N_para+1:end), [6, 6+N_para]);
                    zeros(N_para, 6) eye(N_para, N_para)];
%                 STM = STM .* scale_STM;         % 无量纲化
                H = [dYdX(:, 1:6+N_para) * STM, dYdX(:, 6+N_para+1:end)];
%                 yc = y_c(y_start:i-1) / LEN;    % 无量纲化
                yc = y_c(y_start:i-1);
                Pi = diag(P);
                H0yi{clk_mi, 1} = H;
                H0yi{clk_mi, 2} = yc;
                H0yi{clk_mi, 3} = Pi;
                Hi = 1 * ones(size(dYdX, 1), 1);
%                 B = B + H' * Pi * H - (H' * Pi * Hi) / (Hi' * Pi * Hi) * (Hi' * Pi * H);
%                 L = L + H' * Pi * yc- (H' * Pi * Hi) / (Hi' * Pi * Hi) * (Hi' * Pi * yc);
%                 Hi = eye(size(dYdX, 1));
                B = B + H' * Pi * H - (H' * Pi * Hi) * ((Hi' * Pi * Hi) \ (Hi' * Pi * H));
                L = L + H' * Pi * yc- (H' * Pi * Hi) * ((Hi' * Pi * Hi) \ (Hi' * Pi * yc));
            end
            
            if i == N_all + 1
                break
            end
            
            clk_mi = clk_mi + 1;
            fepo = epoch;
            
            % 叠加矩阵初始化及第一行元素
            y_start = i;
            dYdamb = zeros(1, N_phase_ambiguity);
            if ~isnan(combinedata(i,3))
                dYdamb(amb_list == combinedata(i,3)) = 1;
%                 C = [C; combinedata(i,3)];
            end
            dYdX = [rrhoi(i,:) / ri(i), ...
                0, 0, 0, ...
                zeros(1, N_para), ...
                dYdamb];
            P = combinedata(i, 11);
        else
            dYdamb = zeros(1, N_phase_ambiguity);
            if ~isnan(combinedata(i,3))
                dYdamb(amb_list == combinedata(i,3)) = 1;
%                 C = [C; combinedata(i,3)];
            end
            dYdXnew = [rrhoi(i,:) / ri(i), ...
                0, 0, 0, ...
                zeros(1, N_para), ...
                dYdamb];
            dYdX = [dYdX; dYdXnew];
            P = [P; combinedata(i, 11)];
        end
    end
    if rank(B) < size(B,1)
        B = B + eye(size(B,1));
    end
    DX = B \ L;
    
    % 逐历元计算钟差
    DXclk = zeros(N_clk_error, 1);
    for i = 1: N_clk_error
        N_row = size(H0yi{i, 1}, 1);
        Hii = ones(N_row, 1);
        DXclk(i) = (Hii' * H0yi{i, 3} * Hii) \ (Hii' * H0yi{i, 3} * H0yi{i, 2} - Hii' * H0yi{i, 3} * H0yi{i, 1} * DX);
    end
%     scale_DX = [LEN * ones(3,1); LEN * TAU * ones(3,1); ones(N_para,1); LEN * ones(N_phase_ambiguity,1)];
%     DX = DX .* scale_DX;
%     DXclk = DXclk * LEN;
    toc
    
    %% 收敛判断
    % 计算均方根差
    RMSE = RMSE_new;
    RMSE_new = sqrt(sum(y_c .^2 .* combinedata(:,11)) / N_all);
    
    %% 迭代赋值
    init_state = init_state + [DX(1:6+N_para); DX(6+N_para+1:6+N_para+N_phase_ambiguity); DXclk];
    X_correct = X_correct + DX(1:6+N_para);
    
    if ((abs(RMSE_new - RMSE) <= 0.01) && (abs(RMSE_new - RMSE) >= 0.008)) || (N_LOOP >= 2)
        residual_limit = 1;
    elseif abs(RMSE_new - RMSE) < 0.008
        residual_limit = 2;
    else
        residual_limit = 0;
    end
    
    % 3 sigma 准则，剔除野值
    yc_p = y_c(isnan(combinedata(:,3)));
    yc_l = y_c(~isnan(combinedata(:,3)));
    yc_p_prn = combinedata(isnan(combinedata(:,3)), 4);
    yc_l_prn = combinedata(~isnan(combinedata(:,3)), 4);
    yc_p_ambnum = combinedata(isnan(combinedata(:,3)), 3);
    yc_l_ambnum = combinedata(~isnan(combinedata(:,3)), 3);
    sigma3_p = std(yc_p) * 3;
    sigma3_l = std(yc_l) * 3;
    
    disp(abs(RMSE_new - RMSE));
    disp(['RMS: ',num2str(RMSE_new)]);
    disp(['RMS_L: ',num2str(std(yc_l))]);
    disp(['RMS_P: ',num2str(std(yc_p))]);
%     disp(X_correct);

    t_l = combinedata(~isnan(combinedata(:,3)), 12);
    t_p = combinedata(isnan(combinedata(:,3)), 12);
    dlmwrite(['./output/resid/resid_p_',num2str(N_LOOP),'.dat'], [t_p yc_p yc_p_ambnum])
    dlmwrite(['./output/resid/resid_l_',num2str(N_LOOP),'.dat'], [t_l yc_l yc_l_ambnum])
    dlmwrite(['./output/resid/info_all_',num2str(N_LOOP),'.dat'], [combinedata y_c])
    clklist = zeros(N_clk_error,2);
    clklist(:,2) = init_state(6+N_para+N_phase_ambiguity+1:end);
    for i = 1:N_clk_error
        clkt = combinedata(combinedata(:,2) == clk_num(i),12);
        clklist(i,1) = clkt(1);
    end
    clklist = sortrows(clklist, 1);
    
    cla
    clf
%     subplot(2,2,1)
%     plot(t_p, yc_p, '.')
%     title('伪距残差(m)')
%     xlabel('time(s)')
    title(["第",num2str(N_LOOP),"次迭代"])
    vector_color_gps=rand(32,3);
    for i=1:32
        satid=i;
        indexsat=yc_p_prn==satid;
        tp=t_p(indexsat);
        resp=yc_p(indexsat);
        subplot(2,2,1)
        plot(tp ,resp,'.','Color',vector_color_gps(i,:));
        grid on;hold on;
        ylabel('伪距残差 m')
        
        indexsat=yc_l_prn==satid;
        tp=t_p(indexsat);
        resl=yc_l(indexsat);
        subplot(2,2,3)
        plot(tp ,resl,'.','Color',vector_color_gps(i,:));
        grid on;hold on;
        ylabel('相位残差 m')
    end

%     
%     subplot(2,2,3)
%     plot(t_l, yc_l, '.')
%     title('相位残差(m)')
%     xlabel('time(s)')
    
    subplot(2,2,[2,4])
    plot(clklist(:,1), clklist(:,2)/c*1e9, '.')
    title('钟差序列(ns)')
    xlabel('time(s)')
    
    figname = ['./output/resid/',num2str(N_LOOP), '.fig'];
    saveas(gca, figname)
    figname = ['./output/resid/',num2str(N_LOOP), '.png'];
    saveas(gca, figname)
    drawnow

    if (N_LOOP > 1) && (residual_limit == 1)
        list_num = [];
        i = 1;
        while i <= N_all
            if isnan(combinedata(i,3))
                if abs(y_c(i)) > sigma3_p
                    list_num = [list_num; i];
                end
            else
                if abs(y_c(i)) > sigma3_l
                    list_num = [list_num; i];
                end
            end
            i = i + 1;
        end
        if ~isempty(list_num)
            ode_flag = true;
            combinedata(list_num,:) = [];

            % 剔除卫星数少于4的历元
            epochs = unique(combinedata(:,2));
            N_list_num = size(list_num, 1);
            for i = 1: size(epochs, 1)
                EPO = epochs(i);
                N_EPO = length(find(combinedata(:,2)==EPO));
                if N_EPO < 4
                    combinedata(combinedata(:,2)==EPO, :) = [];
                    N_list_num = N_list_num + N_EPO;
                end
            end
            clear EPO N_EPO;

            fprintf('已剔除%d个野值！\n', N_list_num);
            continue
        else
            ode_flag = true;
        end
    end
end
final_orbit = init_orbit(1:6+N_para) + X_correct;
final_clk = clklist;
final_amb = init_state(6+N_para+1:6+N_para+N_phase_ambiguity);
orbit_list = interp1(time, vector_y, (combinedata(1,10):combinedata(end,10)), 'spline');
orbit_list_itrf = zeros(size(orbit_list, 1), 6);
for i = 1 : size(orbit_list, 1)
     eop = TOOL_geteop(combinedata(1,13)-2400000.5+combinedata(1,14)+(i-1)/86400);
     [r_itrf, v_itrf] = gcrs2itrs(orbit_list(i, 1:3)', orbit_list(i, 4:6)', combinedata(1,13), (i-1)/86400, eop.xp, eop.yp, eop.LOD, eop.UT1_UTC, eop.dX, eop.dY);
     orbit_list_itrf(i, :) = [r_itrf' v_itrf'];
end
sateph2sp3(combinedata(1,13)+18/86400, orbit_list_itrf, 1);
dlmwrite('./output/clk.dat', final_clk);

RMSE = RMSE_new;
y_error = y_c;
end
