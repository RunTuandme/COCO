function [r_GCRS, v_GCRS, HG] = itrs2gcrs(r_ITRF, v_ITRF, JD1, JD2, xp, yp, LOD, UT1_UTC, dX, dY)
% DESCRIPTION:     ITRS to GCRS.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           r_GCRS, GCRS position, unit: m.
%                  v_GCRS, GCRS velocity, unit: m/s.
%                  JD1, JD2, the Julian Date.
%                  xp, yp, the polar motion.
%                  LOD, LOD, the length of day.
%                  UT1_UTC, UT1_UTC, the UT1-UTC.
%                  dX, dY, the nutation in longitude and latitude. position.

        % TT time
    JDTT1 = JD1;
    JDTT2 = JD2 + 69.184/86400;
    
    w_plus = 7.292115146706979e-5 * (1 - LOD / 86400);
    vector_w_plus = [0; 0; w_plus];
    
    % HG: GCRS-to-TIRS
    BPNR_T = sofa_C2TIRS(JD1, JD2, UT1_UTC, dX, dY);
    
    % W: TIRS-to-ITRF
    sp = iauSp00(JDTT1, JDTT2);
    WT = Rx(-yp) * Ry(-xp) * Rz(sp);
    
    r_TIRS = WT' * r_ITRF;

    r_GCRS = BPNR_T' * r_TIRS;    
    v_GCRS = BPNR_T' * (WT' * v_ITRF + cross(vector_w_plus, r_TIRS));
    
    HG = WT * BPNR_T;
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

