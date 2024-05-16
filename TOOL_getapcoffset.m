function [f3d, satfixed2earthfixed] = TOOL_getapcoffset(JD1, JD2, leopos, HG)
% DESCRIPTION:     Get phase center offset of receiver. 
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    global apcoffset 

    JDTT1 = JD1;
    JDTT2 = JD2 + 69.184/86400;
    PVsun = iauEpv00(JDTT1, JDTT2);
    rSun = PVsun(1:3);
    rleo2sun = HG * rSun - leopos;
    ez = -leopos / norm(leopos);
    ey = cross(ez, rleo2sun/norm(rleo2sun));
    ex = cross(ey, ez);
    satfixed2earthfixed = [ex ey ez];
    
    f3d = satfixed2earthfixed * apcoffset;
end