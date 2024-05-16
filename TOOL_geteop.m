function eop = TOOL_geteop(mjd)
% DESCRIPTION:     Get EOP.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    global interxp
    global interyp
    global interUT1_TUC
    global interLOD
    global interdX
    global interdY
    eop.xp = interxp(mjd) / 3600 * pi / 180;
    eop.yp = interyp(mjd) / 3600 * pi / 180;
    eop.UT1_UTC = interUT1_TUC(mjd) / 86400;
    eop.LOD = interLOD(mjd) * 1e-3;
    eop.dX = interdX(mjd) * 1e-3 / 3600 * pi / 180;
    eop.dY = interdY(mjd) * 1e-3 / 3600 * pi / 180;
end