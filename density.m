function [rho] = density(x, y, z, jd, HG)
% DESCRIPTION:     nrlmsise00 density.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           x, y, z, satellite position.
%                  jd, the Julian Date.
%                  HG, rotation matrix from GCRF to ITRF.

    [f107A,f107,ap] = deal(150,150,6);

    ap_a = struct('aa',zeros(7,1));
    ap_a.aa(1) = ap;
    ap_a.aa(2) = ap;
    ap_a.aa(3) = ap;
    ap_a.aa(4) = ap;
    ap_a.aa(5) = ap;
    ap_a.aa(6) = ap;
    ap_a.aa(7) = ap;
    
    [year,month,day,hour,minut,sec] = jd2ce(jd);
%     [latitude,longtitude,h] = starpoint(x,y,z,jd);
%     HG = sofa_C2T(jd);
    r_WGS84 = HG * [x; y; z];
    wx = r_WGS84(1);
    wy = r_WGS84(2);
    wz = r_WGS84(3);
    [latitude,longtitude,h] = pv2ob(wx, wy, wz);
    h = h / 1e3;    % km
    lst= suntime(hour,minut,sec,longtitude);

    input = struct('year',0,'doy',0,'sec',0,'alt',0,'g_lat',0,'g_long',0,'lst',0,'f107A',0,'f107',0,'ap',0,'ap_a',ap_a);
    doy = finddays(year, month, day, hour, minut, sec);
    input.doy = doy;
    input.year = year; % without effect
    input.sec = hour * 3600 + minut * 60 + sec;
    input.alt = h;
%     input.alt = 1109.6339067673;
    input.g_lat = latitude;
%     input.g_lat =25.9542037133;
    input.g_long = longtitude;
%     input.g_long = 96.4205061376;
    input.lst = lst;
%     input.lst = hour + minut / 60 + sec / 3600;
    input.f107A = f107A;
    input.f107 = f107;
    input.ap = ap;
    input.ap_a = ap_a;
    output = nrlmsise00(input);
    rho = output.d(6);      % g/cm^3
    rho = rho * 1e3;        % kg/m^3
end

