
function [final_orbit, JD2_init, orbits] = IOD(orbit)
% DESCRIPTION:     Compute the initial orbit.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% NOTES:           orbit(format), JD1 | JD2 | X | Y | Z | clk | RMS | t0
orbits = orbit(~any(isnan(orbit), 2), :);
orbits = orbits(:,:);

%% 
N_orbits = size(orbits, 1);
t_end = orbits(end-1, end) - orbits(1, end);
tint = orbits(:,end)-orbits(1,end);
%% 地固系转惯性系

for i = 1:N_orbits
    JD1 = orbits(i,1);
    JD2 = orbits(i,2);
    HG = sofa_C2T(JD1, JD2);
    r_GCRS = HG' * orbits(i, 3:5)';
    orbits(i,3:5) = r_GCRS';
end

px = interp1(tint, orbits(:,3), 'spline','pp');
py = interp1(tint, orbits(:,4), 'spline','pp');
pz = interp1(tint, orbits(:,5), 'spline','pp');
velocity = [px.coefs(:,3) py.coefs(:,3) pz.coefs(:,3)];
init_velocity = velocity(1,:);
orbits(end,:) = [];

%% 初始化
init_state = [orbits(1,3:5)'; init_velocity'];               %   初轨计算值
init_orbit = [orbits(1,3:5)'; init_velocity'];              
global JD1_pod
JD1_pod = orbits(1,1);
JD2_init = orbits(1,2);
% t_init = orbits(1,1);                      %   初轨历元时刻

RMSE = 1;
RMSE_new = 0;

CdAM = 0.044;
CrAM = 0.02;

% 大气光压参数分段
global Cd_tinter
global Cr_tinter
global num_Cd
global num_Cr
global exp_switch
Cd_tinter = 3600 * 3;  % second
Cr_tinter = 3600 * 12; % second
% Cd_tinter = inf;  % second
% Cr_tinter = inf; % second
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
H_init = [X_correct; H_init(:)]; % 8状态量 + 6*8 状态转移矩阵

% dYdX = eye(6+N_para);
dYdX = [eye(3) zeros(3,3+N_para)];

init_state = [init_state; CdAM*ones(num_Cd,1); CrAM*ones(num_Cr,1); CC; Cn; Sn];   % 当前估计值
init_orbit = [init_orbit; CdAM*ones(num_Cd,1); CrAM*ones(num_Cr,1); CC; Cn; Sn];   % 先验值

opts = odeset('Reltol',1e-12,'AbsTol',1e-13,'Stats','off','NormControl','on');

N_LOOP = 0;
%% 计算轨道积分&状态转移矩阵
while((abs(RMSE_new - RMSE) >= 0.01) && N_LOOP < 20)
    
%     init_state = data_j2000(1,4:9)';               %   初轨计算值
%     t_init = data_j2000(1,1);                      %   初轨历元时刻
%     N_data_j2000 = size(data_j2000,1);
%     t_end = floor((data_j2000(N_data_j2000, 1) - t_init) * 86400);
    
    N_LOOP = N_LOOP + 1;
    fprintf('开始第%d次迭代...\n', N_LOOP);
    
    tic
    
    H_init(1:length(init_state)) = init_state;

    [time, vector_y] = ode113(@variationalEquation, [0 t_end], H_init, opts, JD2_init);
%     vector_yin = interp1(time, vector_y, orbits(:,end) - orbits(1,end), 'spline');

    t_inter = orbits(:,end) - orbits(1,end);
    N_inter = length(t_inter);
    N_ycol = size(vector_y, 2);
    vector_yin = zeros(N_inter, N_ycol);
    for i = 1:N_ycol
        ilist = interp_larange_list(time, vector_y(:,i), t_inter);
        vector_yin(:,i) = ilist;
    end

    aN = size(vector_yin,1);
    STM = cell(aN, 1);
    for i = 1:aN
        y = vector_yin(i,:)';
        STM{i,1} = [reshape(y(6+N_para+1:end), [6, 6+N_para]);
                    zeros(N_para, 6) eye(N_para, N_para)];
    end
    toc

    %% 计算权矩阵
    W = cell(aN, 1);
    for i = 1:aN
        W{i,1} = eye(3);
    end
    
    %% 计算加权法方程
    NW = cell(aN, 1);
    YW = cell(aN, 1);
    
    % 计算残差
    y_c = orbits(:,3:5) - vector_yin(:,1:3);
  
    % 计算法方程
    for i = 1:aN
        B = dYdX * STM{i};
        NW{i,1} = B' * W{i,1} * B;
        YW{i,1} = B' * W{i,1} * y_c(i,:)';
    end

    %% 批处理
    sumN = 0;
    sumY = 0;
    for i = 1:aN
        sumN = sumN + NW{i,1};
        sumY = sumY + YW{i,1};
    end
    if rank(sumN) < size(sumN, 1)
        sumN = sumN + eye(size(sumN, 1));
    end
    DX = sumN \ sumY; 

    %% 收敛判断
    % 计算均方根差
    RMSE = RMSE_new;
    sum_error = 0;
    for i = 1:aN
        sum_error = sum_error + y_c(i,1:3)  * y_c(i,1:3)';
    end
    RMSE_new = sqrt(sum_error / (aN) / 3);
    
    %% 迭代赋值
    init_state = init_state + DX;
    X_correct = X_correct + DX;
    disp(abs(RMSE_new - RMSE));
    disp(RMSE_new);
    disp(X_correct);
    
    cla
    plot(orbits(:,2), y_c(:,1:3), '.')
    drawnow
    
    % 3 sigma 准则，剔除野值
%     sigma = sqrt(sum(y_c.^2) / aN);
%     sigma3 = sigma * 3;
    sigma3 = std(y_c) * 3;
    
    list_num = [];
    if N_LOOP > 1 || (abs(RMSE_new - RMSE) < 1)
        abs_yc = abs(y_c);
        for i = 1:aN
            if abs_yc(i,1) > sigma3(1) || abs_yc(i,2) > sigma3(2) || abs_yc(i,3) > sigma3(3) 
                list_num = [list_num; i];
            end
        end
        if list_num(1) == 1
            list_num(1) = [];
        end
        if ~isempty(list_num)
%              ti = orbits(list_num,7);
%             for kk = 1:length(ti)
%                 obsi = find(obsdata.epoch(:,3)==ti(kk));
%                 obsdata.P1(obsi,:) = NaN;
%                 obsdata.P2(obsi,:) = NaN;
%                 obsdata.P(obsi,:) = NaN;
%                 obsdata.L1(obsi,:) = NaN;
%                 obsdata.L2(obsi,:) = NaN;
%                 obsdata.L(obsi,:) = NaN;    
%             end
            orbits(list_num,:) = [];
            disp("已剔除IOD数据野值")
            continue
        end
    end
    
end

final_orbit = init_orbit + X_correct;

end