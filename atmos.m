function [drag] = atmos(x, y, z, vx, vy, vz, rho, CdAM)
% DESCRIPTION:     Calculate drag acceleration.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           x, y, z, vx, vy, vz, satellite position and velocity.
%                  rho, the air density, the unit is kg/m^3.
%                  CdAM, drag coefficient, the unit is m^2/kg.
    omega_earth = 7.292115146706979e-5; % rad/s
    delta_vj2000 = cross([0; 0; omega_earth], [x; y; z]);
    vector_vr = [vx; vy; vz] - delta_vj2000;
    vr = norm(vector_vr);
    drag = -rho * CdAM / 2 * vr * vector_vr;
end