 %% ----- 变分方程右函数 ----- %%
function [dHdt] = variationalEquation(t, H, T)
% DESCRIPTION:     Function for numerical integration.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           t, integration time.
%                  H, integration vector.
%                  T, begining time of integration.
% OUTPUT:          dHdt, derivative of H.
% NOTES:
%                  Variational Equation： d2H/dt2 = A H + B dH/dt + C
%                  1. A is the partial derivative of acceleration with respect to position             (3*3)
%                  2. B is the partial derivative of acceleration with respect to velocity             (3*3)
%                  3. C is the partial derivative of acceleration with respect to parameter            (O_3*3, O_3*3, d2r/dβ)
%                  4. H is a vector formed by concatenating the state quantity and the Jacobian matrix (X(6) ,Z(nxn)) 
    
    order = 70;     % 重力场阶数

    % 数据预算
    [x, y, z, vx, vy, vz] = deal(H(1), H(2), H(3), H(4), H(5), H(6));
    global JD1_pod
    global exp_switch
    global Cd_tinter
    global Cr_tinter
    global num_Cd
    global num_Cr
    if num_Cd == 0
        CdAM = 0.044;
    else
        CdAM = H(7 + floor(t / Cd_tinter));
    end
    if num_Cr == 0
        CrAM = 0.02;
    else
        CrAM = H(7 + num_Cd + floor(t / Cr_tinter));
    end
    if exp_switch
        CC = H(7 + num_Cd + num_Cr: 9 + num_Cd + num_Cr);
        Cn = H(10 + num_Cd + num_Cr: 12 + num_Cd + num_Cr);
        Sn = H(13 + num_Cd + num_Cr: 15 + num_Cd + num_Cr);
        N_para = num_Cd + num_Cr + 9;
    else
        CC = zeros(3,1);
        Cn = zeros(3,1);
        Sn = zeros(3,1);
        N_para = num_Cd + num_Cr;
    end

    JD1 = JD1_pod;
    JD2 = T + t / 86400;
    % GCRS-to-TIRS-to-ITRF
    HG = sofa_C2T(JD1, JD2);
    
    rho = density(x, y, z, JD1+JD2, HG);
    Mjd_TT = JD1 - 2400000.5 + JD2 + 69.184/86400;
    Mjd_TDB = Mjday_TDB(Mjd_TT);
    [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
          r_Neptune,r_Pluto,r_Moon,r_Sun,~] = JPL_Eph_DE440(Mjd_TDB);
%     [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
%           r_Neptune,r_Pluto,r_Moon, r_Sun] = JPL_Eph_DE405(jd - 2400000.5);
    JPL = [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
          r_Neptune,r_Pluto,r_Moon,r_Sun];

    %  step1. 轨道预报
    [a, dc, ds, dc0] = acceleration(JD1, JD2, x, y, z, vx, vy, vz, HG, order, JPL, rho, CdAM, CrAM, CC, Cn, Sn);
    [ax, ay, az] = deal(a(1), a(2), a(3)); 

    %  step2. 计算状态转移矩阵
    D = reshape(H(6+N_para+1:end),[6,6+N_para]);
    FUN = D(1:3,:);
    dFUN = D(4:6,:);

    [A, B, C] = derivative(x, y, z, vx, vy, vz, HG, order, dc, ds, dc0, JPL, rho, CdAM, CrAM, max(ceil(t/Cd_tinter),1), max(ceil(t / Cr_tinter),1), JD1, JD2);
    if num_Cr == 0 && num_Cd == 0 && ~exp_switch
        C = zeros(3,6);
    end

    J = A * FUN + B * dFUN + C; % 3 * 17
    K = [dFUN; J];
    dHdt = [vx; vy; vz; ax; ay; az; zeros(N_para,1); K(:)];
end
