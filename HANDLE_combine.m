function obsdata = HANDLE_combine(obsdata, dcbdata)
% DESCRIPTION:     Combine observation data.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    f1 = 1575.42 * 1E6;   % Hz
    f2 = 1227.60 * 1E6;   % Hz
    
    dcb_G = dcbdata.G;
    dcb_G(isnan(dcb_G)) = 0;
    for i = 1: size(dcb_G,2)
        obsdata.P1(:,i) = obsdata.P1(:,i) + dcb_G(1,i);
        obsdata.P2(:,i) = obsdata.P2(:,i) + dcb_G(2,i);
    end
    
    obsdata.L = (f1^2 * obsdata.L1 - f2^2 * obsdata.L2) / (f1^2 - f2^2);
    obsdata.P = (f1^2 * obsdata.P1 - f2^2 * obsdata.P2) / (f1^2 - f2^2);
end