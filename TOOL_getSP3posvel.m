% activeS:    PRN号（GPS） 如[1 3 4 5 17 24]
%             PRN号+GPS卫星个数（BDS）
%             PRN号+GPS卫星个数+BDS卫星个数（Galileo）
% POS    :    ITRS-based


function [POS, VEL] = TOOL_getSP3posvel(rho, t0, JD1, JD2, activeS, HG)
% DESCRIPTION:     Get SP3 position and velocity.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           rho, distance of receiver to satellite.
%                  t0, second.
%                  JD1, the large part of the Julian Date.
%                  JD2, the small part of the Julian Date.
%                  activeS, Active GNSS satellites in search.
%                  HG, rotation matrix from GCRF to ITRF.
    global interpSP3posF
    global interpSP3velF
    global atxdata
    N_GPS = 32;
    c = 2.99792458E+8;
    f1 = 1575.42 * 1E6;             % Hz
    f2 = 1227.60 * 1E6;             % Hz
    N1 = f1^2 / (f1^2 - f2^2);
    N2 = f2^2 / (f1^2 - f2^2);
    
    POS = NaN(length(activeS),3);
    VEL = NaN(length(activeS),3);
    
%     tlist = cell2mat(interpSP3posF.GridVectors);
%     index = find(tlist <= t0, 1, 'last');
%     % 检查找到的结果
%     if isempty(index)
%         error('GNSS插值时间位于序列外\n');
%     else
%         for i = 1:length(activeS)
% 
%             xi = lagrange_interp(tlist(index-5:index+5), interpSP3posF.Values(index-5:index+5, 1, activeS(i)), t0 - rho(i)/c, 11);
%             yi = lagrange_interp(tlist(index-5:index+5), interpSP3posF.Values(index-5:index+5, 2, activeS(i)), t0 - rho(i)/c, 11);
%             zi = lagrange_interp(tlist(index-5:index+5), interpSP3posF.Values(index-5:index+5, 3, activeS(i)), t0 - rho(i)/c, 11);
%             
%             vxi = lagrange_interp(tlist(index-5:index+5), interpSP3velF.Values(index-5:index+5, 1, activeS(i)), t0 - rho(i)/c, 11);
%             vyi = lagrange_interp(tlist(index-5:index+5), interpSP3velF.Values(index-5:index+5, 2, activeS(i)), t0 - rho(i)/c, 11);
%             vzi = lagrange_interp(tlist(index-5:index+5), interpSP3velF.Values(index-5:index+5, 3, activeS(i)), t0 - rho(i)/c, 11);
% 
%             POS(i,:) = [xi, yi, zi];
%             VEL(i,:) = [vxi,vyi,vzi];
%         end
%     end
    
    for i = 1:length(activeS)
        try
        P = interpSP3posF(t0 - rho(i) / c);
        catch
            a=0;
        end
        try
        POS(i,:) = P(:, :, activeS(i));
        catch
            a=0;
        end
        V = interpSP3velF(t0 - rho(i) / c);
        VEL(i,:) = V(:, :, activeS(i));
    end
        
    JDTT1 = JD1;
    JDTT2 = JD2 + 69.184/86400;
    PVsun = iauEpv00(JDTT1, JDTT2);
    rSun = PVsun(1:3);
    
    
    for i = 1:length(activeS)
        if i <= N_GPS
            dcb1 = atxdata(1:3, activeS(i)) / 1e3; % 转换到m
            dcb2 = atxdata(4:6, activeS(i)) / 1e3; % 转换到m
            dcb = N1 * dcb1 - N2 * dcb2;
            
%             rgnss2sun = HG * rSun - POS(i,:)';
%             ez = -POS(i,:)' / norm(POS(i,:));
%             ey = cross(ez, rgnss2sun/norm(rgnss2sun));
%             ex = cross(ey, ez);
%             satfixed2earthfixed = [ex ey ez];

            POS_GCRS = HG' * POS(i,:)';
            rgnss2sun = rSun - POS_GCRS;
            ez = -POS_GCRS / norm(POS_GCRS);
            ey = cross(ez, rgnss2sun/norm(rgnss2sun));
            ex = cross(ey, ez);
            satfixed2inertial = [ex ey ez];

%             sp3modified = satfixed2inertial * dcb;
%             POS(i,:) = POS(i,:) + sp3modified';
            POSmodified = satfixed2inertial * dcb;
            POS_GCRS = POS_GCRS + POSmodified;
            POS(i,:) = (HG * POS_GCRS)';
        end
    end
%     sp3modified = satfixed2earthfixed * dcb;
%     POS = POS - permute(sp3modified, [3 1 2]);
%     POS = POS - sp3modified';
end