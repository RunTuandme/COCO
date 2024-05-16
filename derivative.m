function [A, B, C] = derivative(x, y, z, vx, vy, vz, HG, order, dc, ds, dc0, JPL, rho, CdAM, CrAM, Cdn, Crn, JD1, JD2)
% DESCRIPTION:     Calculate derivative.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           JD1, the large part of the Julian Date.
%                  JD2, the small part of the Julian Date.
%                  HG, rotation matrix from GCRF to ITRF.
%                  x, y, z, vx, vy, vz, satellite position and velocity.
%                  order, order of the gravitational field model.
%                  dc, ds, dc0, coefficients of tide model.
%                  JPL, JPL results, the format is the gather of the output of the function JPL_Eph_DE440.m. e.g. [r_Mercury,r_Venus].
%                  rho, the air density, the unit is kg/m^3.
%                  CdAM, drag coefficient, the unit is m^2/kg.
%                  CrAM, reflectivity coefficient, the unit is m^2/kg.
%                  Cdn, Crn, number of Cd and Cr.
% NOTES:
%                  Variational Equation： d2H/dt2 = A H + B dH/dt + C
%                  1. A is the partial derivative of acceleration with respect to position (3*3)
%                  2. B is the partial derivative of acceleration with respect to velocity (3*3)
%                  3. C is the partial derivative of acceleration with respect to parameter(O_3*3, O_3*3, d2r/dβ)
global GM
global ae
global gravity_type
r = norm([x;y;z]);

%% ------  A  ------ %%
% ----- 中心天体 ----- %
A1 = GM / r^3 * [3*x^2/r^2-1 3*x*y/r^2   3*x*z/r^2;...
                 3*x*y/r^2   3*y^2/r^2-1 3*y*z/r^2;...
                 3*x*z/r^2   3*y*z/r^2   3*z^2/r^2-1];

% ----- 非球形引力 ----- %
global clm0
global clm
global slm

if gravity_type == 3
    global tc0 tc ts acosc0 acosc acoss asinc0 asinc asins
    mjd = JD1 - 2400000.5 + JD2;
    mjd0 = 55197;
    T = 365.25;
    dt = (mjd - mjd0) / T;
    clm0cp = clm0 + tc0 * dt + acosc0 * cos(2*pi*dt) + asinc0 * sin(2*pi*dt) + dc0';
    clmcp = clm + tc * dt + acosc * cos(2*pi*dt) + asinc * sin(2*pi*dt) + dc;
    slmcp = slm + ts * dt + acoss * cos(2*pi*dt) + asins * sin(2*pi*dt) + ds;
else
    clm0cp = clm0 + dc0';
    clmcp = clm + dc;
    slmcp = slm + ds;
end

% ae = 1;
% GM = 1;
XYZ = HG * [x; y; z]; % 转为地固坐标系
X = XYZ(1);
Y = XYZ(2);
Z = XYZ(3);
[lam,  phi, R] = d2l(X, Y, Z);
ae_R = ae / R;
sphi = sin(phi);
cphi = cos(phi);
dRdX = X / R;
dRdY = Y / R;
dRdZ = sphi;
R_ = norm([X;Y]);
dphidX = -dRdX * dRdZ / R_;
dphidY = -dRdY * dRdZ / R_;
dphidZ = R_ / r^2;
dlamdX = -Y / (X^2 + Y^2);
dlamdY = -X / (X^2 + Y^2);
dlamdZ = 0;
% [pl, dpl] = average_Pl_u(sphi, order);
[pl, dpl] = average_Pl_u_mex(sphi, order);
% dpl = dpl * cos(phi);
% [plm, dplm] = average_Plm_u(sphi, order);
[plm, dplm] = average_Plm_u_mex(sphi, order);
% dplm = dplm * cos(phi);
[cosm, sinm] = cos_sin_m(lam, order);
% Tlm
ncosm = repmat(cosm, order, 1);
nsinm = repmat(sinm, order, 1);
Tlm  = (clmcp(1:order, 1:order) .* ncosm + slmcp(1:order, 1:order) .* nsinm) .* ae_R.^meshgrid((1:order),(1:order))';
TClm = (slmcp(1:order, 1:order) .* ncosm - clmcp(1:order, 1:order) .* nsinm) .* ae_R.^meshgrid((1:order),(1:order))';
% d2VdR2
d2VdR2_1 = (2:order+1) .* (3:order+2) .* clm0cp(1:order) .* ae_R.^(1:order) * pl(2:order+1)';
d2VdR2_2 = sum(sum( meshgrid((2:order+1),(1:order))' .* meshgrid((3:order+2),(1:order))' .* plm(1:order,1:order) .* Tlm ));
d2VdR2 = GM * (d2VdR2_1 + d2VdR2_2) / R^3;
% d2VdRdphi
d2VdRdphi_1 = (2:order+1) .* clm0cp(1:order) .* ae_R.^(1:order) * dpl' * cphi;
d2VdRdphi_2 = sum(sum( meshgrid((2:order+1),(1:order))' .* Tlm .* dplm )) * cphi;
d2VdRdphi = -GM * (d2VdRdphi_1 + d2VdRdphi_2) / R^2;
% d2VdRdlam
d2VdRdlam = -GM * sum(sum( meshgrid((1:order),(1:order)) .* meshgrid((2:order+1),(1:order))' .* plm(1:order,1:order) .* TClm )) / R^2;
% d2Vdphidlam
d2Vdphidlam = GM * sum(sum( meshgrid((1:order),(1:order)) .* TClm .* dplm )) * cphi / R;
% d2Vdlam2
d2Vdlam2 = -GM * sum(sum( meshgrid((1:order),(1:order)).^2 .* plm(1:order,1:order) .* Tlm )) / R;
% dVdR
dVdR_1 = (2:order+1) .* clm0cp(1:order) .* ae_R.^(1:order) * pl(2:order+1)';
dVdR_2 = sum(sum( meshgrid((2:order+1),(1:order))' .* plm(1:order,1:order) .* Tlm ));
dVdR = -GM * (dVdR_1 + dVdR_2) / R^2;
% dVdphi
dVdphi_1 = clm0cp(1:order) .* ae_R.^(1:order) * dpl' * cphi;
dVdphi_2 = sum(sum( dplm .* Tlm )) * cphi;
dVdphi = GM * (dVdphi_1 + dVdphi_2) / R;
% dVdlam
dVdlam = GM * sum(sum( plm(1:order,1:order) .* TClm .* meshgrid((1:order),(1:order)) )) / R;
% d2Vdphi2
d2Vdphi2 = -2 * R * dVdR - R^2 * d2VdR2 + tan(phi) * dVdphi - d2Vdlam2 / cos(phi)^2;

% dVdX = dVdR * dRdX + dVdphi * dphidX + dVdlam * dlamdX;
% dVdY = dVdR * dRdY + dVdphi * dphidY + dVdlam * dlamdY;
% dVdZ = dVdR * dRdZ + dVdphi * dphidZ;
ddVdXdR = dRdX * d2VdR2 + dphidX * (d2VdRdphi - dVdphi/R) + dlamdX * (d2VdRdlam - dVdlam/R);
ddVdXdphi = dRdX * (d2VdRdphi - dVdphi/R) + dphidX * (d2Vdphi2 + R*dVdR) + dlamdX * (d2Vdphidlam + tan(phi)*dVdlam);
ddVdXdlam = dRdX * d2VdRdlam + dphidX * d2Vdphidlam + dlamdX * d2Vdlam2 - dRdY * dVdR - dphidY * dVdphi - dlamdY * dVdlam;
ddVdYdR = dRdY * d2VdR2 + dphidY * (d2VdRdphi - dVdphi/R) + dlamdY * (d2VdRdlam - dVdlam/R);
ddVdYdphi = dRdY * (d2VdRdphi - dVdphi/R) + dphidY * (d2Vdphi2 + R*dVdR) + dlamdY * (d2Vdphidlam + tan(phi)*dVdlam);
ddVdYdlam = dRdY * d2VdRdlam + dphidY * d2Vdphidlam + dlamdY * d2Vdlam2 + dVdR * dRdX + dVdphi * dphidX + dVdlam * dlamdX;
ddVdZdR = dRdZ * d2VdR2 + dphidZ * (d2VdRdphi - dVdphi/R);
ddVdZdphi = dRdZ * (d2VdRdphi - dVdphi/R) + dphidZ * (d2Vdphi2 + R * dVdR);
ddVdZdlam = dRdZ * d2VdRdlam + dphidZ * d2Vdphidlam;
d2VdXdX = ddVdXdR * dRdX + ddVdXdphi * dphidX + ddVdXdlam * dlamdX;
d2VdXdY = ddVdXdR * dRdY + ddVdXdphi * dphidY + ddVdXdlam * dlamdY;
d2VdXdZ = ddVdZdR * dRdX + ddVdZdphi * dphidX + ddVdZdlam * dlamdX;
d2VdYdY = ddVdYdR * dRdY + ddVdYdphi * dphidY + ddVdYdlam * dlamdY;
d2VdYdZ = ddVdZdR * dRdY + ddVdZdphi * dphidY + ddVdZdlam * dlamdY;
d2VdZdZ = ddVdZdR * dRdZ + ddVdZdphi * dphidZ + ddVdZdlam * dlamdZ;
a1 = d2VdXdX;
a2 = d2VdXdY;
a3 = d2VdXdZ;
a4 = d2VdXdY;
a5 = d2VdYdY;
a6 = d2VdYdZ;
a7 = d2VdXdZ;
a8 = d2VdYdZ;
a9 = d2VdZdZ;
A2 = [a1 a2 a3;a4 a5 a6;a7 a8 a9];
A2 = HG' * A2 * HG;

% ----- 三体引力 ----- %
    r_Mercury = JPL(:,1);
    r_Venus = JPL(:,2);
    r_Mars = JPL(:,4);
    r_Jupiter = JPL(:,5);
    r_Saturn = JPL(:,6);
    r_Uranus = JPL(:,7);
    r_Neptune = JPL(:,8);
    r_Pluto = JPL(:,9);
    r_Moon = JPL(:,10);
    r_Sun = JPL(:,11);
    
    GM_Earth   = 398600.4415e9;         			   % [m^3/s^2]; GGM03C & GGM03S
    GM_Sun     = 132712440041.279419e9; 			   % [m^3/s^2]; DE440
    GM_Moon    = GM_Earth/81.3005682214972154;         % [m^3/s^2]; DE440
    GM_Mercury = 22031.868551e9; 		  			   % [m^3/s^2]; DE440
    GM_Venus   = 324858.592000e9;       			   % [m^3/s^2]; DE440
    GM_Mars    = 42828.375816e9;	      			   % [m^3/s^2]; DE440
    GM_Jupiter = 126712764.100000e9;    			   % [m^3/s^2]; DE440
    GM_Saturn  = 37940584.841800e9;     			   % [m^3/s^2]; DE440
    GM_Uranus  = 5794556.400000e9;      			   % [m^3/s^2]; DE440
    GM_Neptune = 6836527.100580e9;      			   % [m^3/s^2]; DE440
    GM_Pluto   = 975.500000e9;	      			       % [m^3/s^2]; DE440

A3 = dN3tosat(XYZ, r_Sun, GM_Sun) + ...
     dN3tosat(XYZ, r_Moon, GM_Moon) + ...
     dN3tosat(XYZ, r_Mercury, GM_Mercury) + ...
     dN3tosat(XYZ, r_Venus, GM_Venus) + ...
     dN3tosat(XYZ, r_Mars, GM_Mars) + ...
     dN3tosat(XYZ, r_Jupiter, GM_Jupiter) + ...
     dN3tosat(XYZ, r_Saturn, GM_Saturn) + ...
     dN3tosat(XYZ, r_Uranus, GM_Uranus) + ...
     dN3tosat(XYZ, r_Neptune, GM_Neptune) + ...
     dN3tosat(XYZ, r_Pluto, GM_Pluto);
%     A3 = dN3tosat(XYZ, r_Sun, GM_Sun) + ...
%      dN3tosat(XYZ, r_Moon, GM_Moon);

% ----- 大气阻力 ----- %
% global px py pz
% global mux muy muz
omega_earth = 7.292115146706979e-5; % rad/s
delta_vj2000 = cross([0; 0; omega_earth], [x; y; z]);
vector_vr = [vx; vy; vz] - delta_vj2000;
vrx = vector_vr(1);
vry = vector_vr(2);
vrz = vector_vr(3);
vr = norm(vector_vr);
dadrho = -CdAM / 2 * vr * vector_vr;
% % dadvr
% x_scale = (x-mux(1))/mux(2);
% y_scale = (y-muy(1))/muy(2);
% z_scale = (z-muz(1))/muz(2);
% drhodx = px * [6*x_scale^5; 5*x_scale^4; 4*x_scale^3; 3*x_scale^2; 2*x_scale; 1; 0] / mux(2);
% drhody = py * [6*y_scale^5; 5*y_scale^4; 4*y_scale^3; 3*y_scale^2; 2*y_scale; 1; 0] / muy(2);
% drhodz = pz * [6*z_scale^5; 5*z_scale^4; 4*z_scale^3; 3*z_scale^2; 2*z_scale; 1; 0] / muz(2);
% drhodr = [drhodx drhody drhodz];
drhodr = [0 0 0];
% % dadvr
daxdvrx = -rho * CdAM / 2 * (vr + vrx^2 / vr);
daxdvry = -rho * CdAM / 2 * vrx * vry / vr;
daxdvrz = -rho * CdAM / 2 * vrx * vrz / vr;
daydvrx = daxdvry;
daydvry = -rho * CdAM / 2 * (vr + vry^2 / vr);
daydvrz = -rho * CdAM / 2 * vry * vrz / vr;
dazdvrx = daxdvrz;
dazdvry = daydvrz;
dazdvrz = -rho * CdAM / 2 * (vr + vrz^2 / vr);
dadvr = [daxdvrx daxdvry daxdvrz; ...
         daydvrx daydvry daydvrz; ...
         dazdvrx dazdvry dazdvrz];
dvrdr = [0 omega_earth 0; -omega_earth 0 0; 0 0 0];
A5 = dadrho * drhodr + dadvr * dvrdr;

%
% ----- 光压阻力 ----- %
p0 = 4.5605e-6;         % N/m^2
AU = 1.495978707e11;    % m
vector_r = [x; y; z];
vector_Rs = r_Sun - vector_r;
vector_Rx = vector_Rs(1);
vector_Ry = vector_Rs(2);
vector_Rz = vector_Rs(3);
Rs = norm(vector_Rs);
r = norm(vector_r);

cosPhi = dot(vector_r/r, vector_Rs/Rs);
sinPhi = sqrt(1 - cosPhi^2);

dRsdr = [3*vector_Rx^2/Rs^2-1,       3*vector_Rx*vector_Ry/Rs^2, 3*vector_Rx*vector_Rz/Rs^2;
         3*vector_Rx*vector_Ry/Rs^2, 3*vector_Ry^2/Rs^2-1,       3*vector_Ry*vector_Rz/Rs^2;
         3*vector_Rx*vector_Rz/Rs^2, 3*vector_Ry*vector_Rz/Rs^2, 3*vector_Rz^2/Rs^2-1];

if ( cosPhi <= 0 ) && ( sinPhi < (ae/r) )
    A6 = 0;
else
    A6 = - p0 * (AU / Rs)^2 * CrAM * dRsdr / Rs^3;
end


% A = A1 + A2 + A3 + A4 + A5;
A = A1 + A2 + A3 + A5 + A6;

%% ------  B  ------ %%
B = zeros(3, 3);
% B = dadvr;

%% ------  C  ------ %%
% C = zeros(3, 6);
% 大气弹道系数修正
omega_earth = 7.292115146706979e-5;% rad/s
delta_vj2000 = cross([0; 0; omega_earth], [x; y; z]);
vector_vr = [vx; vy; vz] - delta_vj2000;
vr = norm(vector_vr);
dCddx = -rho / 2 * vr * vector_vr;
global num_Cd
C1 = zeros(3, num_Cd);
C1(:, Cdn) = dCddx;

% 光压参数修正
p0 = 4.5605e-6;         % N/m^2
AU = 1.495978707e11;    % m
vector_r = [x; y; z];
vector_Rs = r_Sun - vector_r;

Rs = norm(vector_Rs);

delta_r = vector_Rs / Rs;
cosPhi = dot(vector_r/r, vector_Rs/Rs);
sinPhi = sqrt(1 - cosPhi^2);
if ( cosPhi <= 0 ) && ( sinPhi < (ae/r) )
    dCrdx = [0; 0; 0];
else
    dCrdx = - p0 * (AU / Rs)^2 * delta_r;
end
global num_Cr
C2 = zeros(3, num_Cr);
C2(:, Crn) = dCrdx;

global exp_switch
if exp_switch
    % 经验力参数修正
    u = 3.9860044150E+14;
    vector_r = [x; y; z];
    vector_v = [vx; vy; vz];
    r = norm(vector_r);
    v = norm(vector_v);
    a = 2 / r - v^2 / u;
    a = 1 / a;
    ecosE = 1 - r / a;
    esinE =(x*vx + y*vy + z*vz)/sqrt(u*a);
    e = sqrt(ecosE^2 + esinE^2);
    cosf = (a-r) / (e*r) - (a*e) / r;
    sinf = a * sqrt(1-e^2) * esinE / (e * r); 
    % 法向量
    vector_n = cross(vector_r, vector_v);
    % RTN 基向量
    er = vector_r / r;
    et = vector_v / v;
    en = vector_n / norm(vector_n);
    % 投影矩阵(J2000-RTN)
    B2RTN = [er'; et'; en'];
    dadCC = B2RTN' * 1e-9;
    dadCn = B2RTN' * cosf * 1e-9;
    dadSn = B2RTN' * sinf * 1e-9;
    C3 = [dadCC dadCn dadSn];
    % C3 = dadCC;

    C = [zeros(3,6) C1 C2 C3];
else
    C = [zeros(3,6) C1 C2];
end

end

%% ----- support function ----- %%
function A = dN3tosat(r_sat, r_M, GM)
    delta_M = r_sat - r_M;
    nr_M = norm(delta_M);
    A = GM/nr_M^3 * [ 3*delta_M(1)^2/nr_M^2-1          3*delta_M(1)*delta_M(2)/nr_M^2   3*delta_M(1)*delta_M(3)/nr_M^2;...
                      3*delta_M(1)*delta_M(2)/nr_M^2   3*delta_M(2)^2/nr_M^2-1          3*delta_M(2)*delta_M(3)/nr_M^2;...
                      3*delta_M(1)*delta_M(3)/nr_M^2   3*delta_M(2)*delta_M(3)/nr_M^2   3*delta_M(3)^2/nr_M^2-1];
end

% % 计算pl,dpl
% function [pl, dpl] = average_Pl_u(u, order)
%     pl = zeros(1, order+1);
%     dpl= zeros(1, order);
%     pl(1) = 1;
%     pl(2) = 1.732050807568877 * u;
%     for l = 2:order
%         part1 = sqrt( (2*l+1)/(2*l-1) );
%         part2 = (2-1/l) * u * pl(l);
%         part3 = sqrt( (2*l-1)/(2*l-3) );
%         part4 = (1-1/l) * pl(l-1);
%         pl(l+1) = part1 * (part2 - part3 * part4);
%     end
%     dpl(1) = 1.732050807568877;
%     for l = 2:order
%         pp1 = l / (1-u^2);
%         pp2 = sqrt( (2*l+1)/(2*l-1) ) * pl(l);
%         pp3 = u * pl(l+1);
%         dpl(l) = pp1 * (pp2 - pp3);
%     end
% end

% % 计算plm，dplm
% function [p, p_] = average_Plm_u(u, order)
%     p = zeros(order+1,order+1);
%     p_ = zeros(order,order);
%     
%     p(1,1) = 1.732050807568877 * sqrt(1-u^2);
%     for l = 2:order+1
%         p(l,l) = sqrt( (2*l+1)/(2*l) ) * sqrt(1-u^2) * p(l-1,l-1);
%     end
%     for l=2:order+1
%         p(l,l-1) = sqrt((2*l+1)) * u * p(l-1,l-1);
%     end
%     for l = 3:order+1
%         for m = 1:l-2
%             p(l,m) = sqrt((2*l+1)*(2*l-1)/(l+m)/(l-m)) * u * p(l-1,m) -...
%                 sqrt((2*l+1)*(l-1+m)*(l-1-m)/(2*l-3)/(l+m)/(l-m)) * p(l-2,m);
%         end
%     end   
%     
%     for l = 1:order
%         for m = 1:l
%             p_(l,m) = 1 / sqrt(1-u^2)...
%                 * (sqrt( (l+m+1)*(l-m) )...
%                 * p(l,m+1) - m * u / sqrt(1-u^2) * p(l,m));
%         end
%     end
% end
