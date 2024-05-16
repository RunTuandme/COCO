function dcbdata = HANDLE_readosb(filename)
% DESCRIPTION:     Read the OSB file.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    global obscodes
    c = 2.99792458E+8;
    fid = fopen(filename, 'r');
    G1 = obscodes.OG{1};
    G2 = obscodes.OG{3};
    C1 = obscodes.OC{1};
    C2 = obscodes.OC{3};
    E1 = obscodes.OE{1};
    E2 = obscodes.OE{3};
    G1 = obscodes.G1;
    G2 = obscodes.G2;
    
    GPSsys = false;
    BDSsys = false;
    GALsys = false;
    if contains(obscodes.sys, 'G')
        GPSsys = true;
        Gnum = obscodes.Gnum;
        dcbdata.G = NaN(2, Gnum);
        gi1 = 0;
        gi2 = 0;
    end
    if contains(obscodes.sys, 'C')
        BDSsys = true;
        Cnum = obscodes.Cnum;
        dcbdata.C = NaN(2, Cnum);
        ci1 = 0;
        ci2 = 0;
    end
    if contains(obscodes.sys, 'E')
        GALsys = true;
        Enum = obscodes.Enum;
        dcbdata.E = NaN(2, Enum);
        ei1 = 0;
        ei2 = 0;
    end
    
    while ~feof(fid)
        str = fgets(fid);
        if contains(str, '-BIAS/SOLUTION')
            break
        end
        if contains(str, 'OSB  G') && GPSsys
            if contains(str(26:28), G1)
                prn = str2double(str(13:14));
                X1bias = str2double(str(83:91)) * c * 1e-9;
                dcbdata.G(1,prn) = X1bias;
                gi1 = gi1 + 1;
            end
            if contains(str(26:28), G2)
                prn = str2double(str(13:14));
                X2bias = str2double(str(83:91)) * c * 1e-9;
                dcbdata.G(2,prn) = X2bias;
                gi2 = gi2 + 1;
            end
            if gi1 == Gnum && gi2 == 0
                GPSsys = false;
                disp(['warning: Lack of ', G2, ' in OSB file'])
            elseif gi2 == Gnum && gi1 == 0
                GPSsys = false;
                disp(['warning: Lack of ', G1, ' in OSB file'])
            elseif gi2 == Gnum && gi1 == Gnum
                GPSsys = false;
            end
        end
        if contains(str, 'OSB  C') && BDSsys
            if contains(str(26:28), C1)
                prn = str2double(str(13:14));
                X1bias = str2double(str(83:91)) * c * 1e-9;
                dcbdata.C(1,prn) = X1bias;
                ci1 = ci1 + 1;
            end
            if contains(str(26:28), C2)
                prn = str2double(str(13:14));
                X2bias = str2double(str(83:91)) * c * 1e-9;
                dcbdata.C(2,prn) = X2bias;
                ci2 = ci2 + 1;
            end
            if ci1 == Cnum && ci2 == 0
                BDSsys = false;
                disp(['warning: Lack of ', C2, ' in OSB file'])
            elseif ci2 == Cnum && ci1 == 0
                BDSsys = false;
                disp(['warning: Lack of ', C1, ' in OSB file'])
            elseif ci2 == Cnum && ci1 == Cnum
                BDSsys = false;
            end
        end
        if contains(str, 'OSB  E') && GALsys
            if contains(str(26:28), E1)
                prn = str2double(str(13:14));
                X1bias = str2double(str(83:91)) * c * 1e-9;
                dcbdata.E(1,prn) = X1bias;
                ei1 = ei1 + 1;
            end
            if contains(str(26:28), E2)
                prn = str2double(str(13:14));
                X2bias = str2double(str(83:91)) * c * 1e-9;
                dcbdata.E(2,prn) = X2bias;
                ei2 = ei2 + 1;
            end
            if ei1 == Enum && ei2 == 0
                GALsys = false;
                disp(['warning: Lack of ', E2, ' in OSB file'])
            elseif ei2 == Enum && ei1 == 0
                GALsys = false;
                disp(['warning: Lack of ', E1, ' in OSB file'])
            elseif ei2 == Enum && ei1 == Enum
                GALsys = false;
            end
        end
        if ~GPSsys && ~BDSsys && ~GALsys
            break
        end
    end
end