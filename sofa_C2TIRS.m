function BPNR_T = sofa_C2TIRS(JD1, JD2, UT1_UTC, dX, dY)
% DESCRIPTION:     Function of SOFA. 
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% NOTES:           GCRS-to-TIRS

    % gain from IERS, or set zero.
    % see: https://www.iers.org/IERS/EN/DataProducts/tools/eop_of_today/eop_of_today_tool.html
%     xp = 0;         % radians
%     yp = 0;         % radians
%     UT1_UTC = 0;
%     global EOPmat
%     mjd_day = floor(jd - 2400000.5);
%     UT1_UTC = EOPmat(EOPmat(:,1)==mjd_day,4) / 86400;
    UT11 = JD1;
    UT12 = JD2 + UT1_UTC;
    JDTT1 = JD1;
    JDTT2 = JD2 + 69.184/86400;
    
    [dpsi, deps] = iauNut00a(JDTT1, JDTT2);
%     [dpsi, deps] = iauNut00a_mex(JDTT1, JDTT2);
%     [M1, M2, epsa] = iauPn06(JDTT, 0, dpsi, deps);
    [epsa, ~, ~, rbp, ~, rbpn] = iauPn00(JDTT1, JDTT2, dpsi, deps);
    
    % /* Transform dX,dY corrections from GCRS to mean of date. */
    v2 = rbpn * [dX; dY; 0];
    ddpsi = v2(1) / sin(epsa);
    ddeps = v2(2);
    dpsi = dpsi + ddpsi;
    deps = deps + ddeps;
    
    N = Rx(-(epsa + deps)) * Rz(-dpsi) * Rx(epsa);
    NPB = N * rbp;
    
%     GST = iauGst06(UT1, 0, JDTT, 0.0, rbpn);
    GST = iauAnp ( iauGmst00 ( UT11, UT12, JDTT1, JDTT2 ) + iauEe00 ( JDTT1, JDTT2, epsa, dpsi ) );
    R = Rz(GST);    
    
    BPNR_T =  R * NPB;
end

function [NPB, PB, eps] = iauPn06(date1, date2, dpsi, deps)
    % Bias-precession Fukushima-Williams angles of date.
    [gamb, phib, psib, eps] = iauPfw06(date1, date2);
    
    % BPN = Rz(gamb) * Rx(phib) * Rz(-(psib + dpsi)) * Rx(-(eps + deps));
    NPB = Rx(-(eps + deps)) * Rz(-(psib + dpsi)) * Rx(phib) * Rz(gamb);
    
    PB = Rx(-eps) * Rz(-psib) * Rx(phib) * Rz(gamb);
end

function [gamb, phib, psib, eps] = iauPfw06(date1, date2)
    DJ00 = 2451545.0;
    DJC = 36525.0;
    
    % Arcseconds to radians
    DAS2R = 4.848136811095359935899141e-6;
    
    % Interval between fundamental date J2000.0 and given date (JC). 
    t = ((date1 - DJ00) + date2) / DJC;
    
    % P03 bias+precession angles.
    gamb = (    -0.052928     + ...
           (    10.556378     + ...
           (     0.4932044    + ...
           (    -0.00031238   + ...
           (    -0.000002788  + ...
           (     0.0000000260 ) ...
           * t) * t) * t) * t) * t) * DAS2R;
    phib = ( 84381.412819     + ...
           (   -46.811016     + ...
           (     0.0511268    + ...
           (     0.00053289   + ...
           (    -0.000000440  + ...
           (    -0.0000000176 ) ...
           * t) * t) * t) * t) * t) * DAS2R;
    psib = (    -0.041775     + ...
           (  5038.481484     + ...
           (     1.5584175    + ...
           (    -0.00018522   + ...
           (    -0.000026452  + ...
           (    -0.0000000148 ) ...
           * t) * t) * t) * t) * t) * DAS2R;
    eps = (84381.406     + ...
          (-46.836769    + ...
          ( -0.0001831   + ...
          (  0.00200340  + ...
          ( -0.000000576 + ...
          ( -0.0000000434) * t) * t) * t) * t) * t) * DAS2R; 
end

function GST = iauGst06(uta, utb, tta, ttb, NPB)
    % Extract CIP coordinates.
    x = NPB(3,1);
    y = NPB(3,2);
    
    % The CIO locator, s.
    s = iauS06(tta, ttb, x, y);
    
    % Greenwich apparent sidereal time.
    era = iauEra00(uta, utb);
    eors = iauEors(NPB, s);
    gst_temp = era - eors;
    GST = rem(gst_temp, 2*pi);
    if (GST < 0) 
       GST = GST + 2*pi;
    end
end

function eo = iauEors(NPB, s)
    % Evaluate Wallace & Capitaine (2006) expression (16). */
    x = NPB(3,1);
    ax = x / (1.0 + NPB(3,3));
    xs = 1.0 - ax * x;
    ys = -ax * NPB(3,2);
    zs = -x;
    p = NPB(1,1) * xs + NPB(1,2) * ys + NPB(1,3) * zs;
    q = NPB(2,1) * xs + NPB(2,2) * ys + NPB(2,3) * zs;
    if (p ~= 0) || (q ~= 0)
        eo = s - atan2(q, p);
    else
        eo = s;
    end
end

function sp = iauSp00(date1, date2)
    DJ00 = 2451545.0;
    DJC = 36525.0;
    
    % Arcseconds to radians
    DAS2R = 4.848136811095359935899141e-6;
    
    % Interval between fundamental date J2000.0 and given date (JC). 
    t = ((date1 - DJ00) + date2) / DJC;
    
    % Approximate s'.
    sp = -47e-6 * t * DAS2R;
end


function C = Rx(theta)
C = [   1              0               0;
        0          cos(theta)      sin(theta);
        0          -sin(theta)     cos(theta)];
end

function C = Ry(theta)
C = [cos(theta)         0           -sin(theta);
        0               1               0;
     sin(theta)         0           cos(theta)];
end

function C = Rz(theta)
C = [cos(theta)     sin(theta)          0;
     -sin(theta)    cos(theta)          0;
        0               0               1];
end