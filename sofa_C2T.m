function HG = sofa_C2T(JD1, JD2)
% DESCRIPTION:     Function of SOFA. 
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% NOTES:           GCRS-to-TIRS-to-ITRF
    MJD1 = JD1 - 2400000.5;
    MJD = MJD1 + JD2;
    eop = TOOL_geteop(MJD);
    
    xp = eop.xp;         % radians
    yp = eop.yp;         % radians
    UT1_UTC = eop.UT1_UTC;  
    dX = eop.dX;
    dY = eop.dY;
    
    JDTT1 = JD1;
    JDTT2 = JD2 + 69.184/86400;
    
    % BPNR_T: GCRS-to-TIRS
    BPNR_T = sofa_C2TIRS(JD1, JD2, UT1_UTC, dX, dY);
    
    % WT: TIRS-to-ITRF
    sp = iauSp00(JDTT1, JDTT2);
    WT = iauRx(-yp) * iauRy(-xp) * iauRz(sp);
    
    HG = WT * BPNR_T;
end