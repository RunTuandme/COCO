function [dfdt, dc, ds, dc0] = acceleration(JD1, JD2, x, y, z, vx, vy, vz, HG, order, JPL, rho, CdAM, CrAM, CC, Cn, Sn)
% DESCRIPTION:     Calculate the acceleration of a satellite.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           JD1, the large part of the Julian Date.
%                  JD2, the small part of the Julian Date.
%                  x, y, z, vx, vy, vz, satellite position and velocity.
%                  order, order of the gravitational field model.
%                  JPL, JPL results, the format is the gather of the output of the function JPL_Eph_DE440.m. e.g. [r_Mercury,r_Venus].
%                  rho, the air density, the unit is kg/m^3.
%                  CdAM, drag coefficient, the unit is m^2/kg.
%                  CrAM, reflectivity coefficient, the unit is m^2/kg.
%                  CC, Cn, Sn, coefficients of the empirical force model.
% Output:          dfdt, All acceleration.
%                  dc, ds, dc0, coefficients of tide model.
    r = [x; y; z];
    nr = norm(r);
    u = 3.9860044150E+14;
    ax = - u * x / nr^3;
    ay = - u * y / nr^3;
    az = - u * z / nr^3;
    g_force_main = [ax;ay;az];
    r_Moon = JPL(:,10);
    r_Sun = JPL(:,11);

    [dc, ds, dc0] = tide(r_Sun, r_Moon, JD1, JD2, HG);
    rb = HG * r;
    [lam,  phi, ~] = d2l(rb(1), rb(2), rb(3));
    [al, alm] = gravityfield(phi, lam, rb(1), rb(2), rb(3), order, dc, ds, dc0, JD1, JD2);
    AL = HG' * al;
    ALM = HG' * alm;

    g_force_J70 = AL + ALM;
    g_force_3body = thirdbody(x, y, z, JPL);
    g_force_light_pressure = light_pressure(x, y, z, r_Sun, CrAM);
    g_force_atmos = atmos(x, y, z, vx, vy, vz, rho, CdAM);
    g_force_rel = relativity(x, y, z, vx, vy, vz);
    g_force_exp = experience(x, y, z, vx, vy, vz, CC, Cn, Sn);
    dfdt = g_force_main + g_force_J70 + g_force_3body + g_force_light_pressure + g_force_atmos + g_force_rel + g_force_exp;
end



