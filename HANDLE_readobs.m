function obsdata = HANDLE_readobs(filename)
% DESCRIPTION:     Read the OBS file.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    global obscodes
    OG = obscodes.OG;
    OC = obscodes.OC;
    OE = obscodes.OE;
    OG_on = false;
    OC_on = false;
    OE_on = false;
    if contains(obscodes.sys, 'G')
        OG_on = true;
    end
    if contains(obscodes.sys, 'C')
        OC_on = true;
    end
    if contains(obscodes.sys, 'E')
        OE_on = true;
    end
    fid = fopen(filename, 'r');
    str = fgets(fid);
    Version = str2double(str(5:9));
    N_eops = 200000;
    obsdata.epoch = NaN(N_eops, 2); % 4包括P1、P2、L1、L2四个观测量
    obsdata.P1 = NaN(N_eops, obscodes.num);
    obsdata.P2 = NaN(N_eops, obscodes.num);
    obsdata.L1 = NaN(N_eops, obscodes.num);
    obsdata.L2 = NaN(N_eops, obscodes.num);
    %%
    if abs(Version-2.20) < eps
        % RINEX 2.20
        i = 0;
        flag_line = false;
        while ~feof(fid)
            str = fgets(fid);
            if contains(str, 'END OF HEADER')
                flag_line = true;
                break;
            end
        end
        while ~feof(fid)
            str = fgets(fid);
            if flag_line
                i = i + 1;
                year = str2double(str(2:3)) + 2000;
                month = str2double(str(5:6));
                day = str2double(str(8:9));
                hour = str2double(str(11:12));
                min = str2double(str(14:15));
                sec = str2double(str(17:26));
                [JD1, JD2] = TOOL_countjd(year ,month, day, hour, min, sec);
                MJD = JD1 - 2400000.5 + JD2;
                MJD_day = floor(MJD);
                sod = round((MJD - MJD_day) * 86400);
                satnum = str2double(str(31:32));
                PRN_str = str(34:end);
                obsdata.epoch(i, :) = [MJD_day sod];
                for sati = 1:satnum
                    PRN = str2double(PRN_str((sati*3-2):(sati*3-1)));
                    str = fgets(fid);
                    obsdata.P1(i, PRN) = str2double(str(50:64));
                    obsdata.P2(i, PRN) = str2double(str(66:80));
                    obsdata.L1(i, PRN) = str2double(str(2:16));
                    obsdata.L2(i, PRN) = str2double(str(18:32));
                    str = fgets(fid);
                end
            end
        end
        % 寻找并去除全为NaN值的行
        rowIsAllNaN = all(isnan(obsdata.P1), 2) & all(isnan(obsdata.P2), 2) & all(isnan(obsdata.L1), 2) & all(isnan(obsdata.L2), 2);
        obsdata.epoch(rowIsAllNaN, :) = [];
        obsdata.P1(rowIsAllNaN, :) = [];
        obsdata.P2(rowIsAllNaN, :) = [];
        obsdata.L1(rowIsAllNaN, :) = [];
        obsdata.L2(rowIsAllNaN, :) = [];
    
    %% 
    elseif abs(Version-3.04) < eps || abs(Version-3.05) < eps
        % RINEX 3.40
        obsdata.OG = zeros(1,4);
        obsdata.OC = zeros(1,4);
        
        if ~OG_on && ~OC_on && OE_on
            error("RINEX3.04读取：未选择GNSS系统");
        end
        i = 0;
        flag_line = false;
        while ~feof(fid)
            str = fgets(fid);
            if contains(str, 'OBS TYPES')
                % ------ GPS ------ %
                if strcmp(str(1), 'G') && OG_on
                    frq_num = str2double(str(5:6));
                    for frqi = 1:frq_num
                        if strcmp(str(4+4*frqi:6+4*frqi), OG{1})
                            obsdata.OG(1) = frqi; % P1
                        elseif strcmp(str(4+4*frqi:6+4*frqi), OG{2})
                            obsdata.OG(2) = frqi; % L1
                        elseif strcmp(str(4+4*frqi:6+4*frqi), OG{3})
                            obsdata.OG(3) = frqi; % P2
                        elseif strcmp(str(4+4*frqi:6+4*frqi), OG{4})
                            obsdata.OG(4) = frqi; % L2
                        else
                            continue
                        end
                    end
                end
                % ------ BDS ------ %
                if strcmp(str(1), 'C') && OC_on
                    frq_num = str2double(str(5:6));
                    for frqi = 1:frq_num
                        if strcmp(str(4+4*frqi:6+4*frqi), OC{1})
                            obsdata.OC(1) = frqi; % P1
                        elseif strcmp(str(4+4*frqi:6+4*frqi), OC{2})
                            obsdata.OC(2) = frqi; % L1
                        elseif strcmp(str(4+4*frqi:6+4*frqi), OC{3})
                            obsdata.OC(3) = frqi; % P2
                        elseif strcmp(str(4+4*frqi:6+4*frqi), OC{4})
                            obsdata.OC(4) = frqi; % L2
                        else
                            continue
                        end
                    end
                end
                % ------ GAL ------ %
                if strcmp(str(1), 'E') && OE_on
                    frq_num = str2double(str(5:6));
                    for frqi = 1:frq_num
                        if strcmp(str(4+4*frqi:6+4*frqi), OE{1})
                            obsdata.OE(1) = frqi; % P1
                        elseif strcmp(str(4+4*frqi:6+4*frqi), OE{2})
                            obsdata.OE(2) = frqi; % L1
                        elseif strcmp(str(4+4*frqi:6+4*frqi), OE{3})
                            obsdata.OE(3) = frqi; % P2
                        elseif strcmp(str(4+4*frqi:6+4*frqi), OE{4})
                            obsdata.OE(4) = frqi; % L2
                        else
                            continue
                        end
                    end
                end
            end
            if contains(str, 'END OF HEADER')
                flag_line = true;
                break;
            end
        end
        
        while ~feof(fid)
            str = fgets(fid);
            if flag_line && strcmp(str(1),'>')
                i = i + 1;
                year = str2double(str(3:6));
                month = str2double(str(8:9));
                day = str2double(str(11:12));
                hour = str2double(str(14:15));
                min = str2double(str(17:18));
                sec = str2double(str(20:29));
                [JD1, JD2] = TOOL_countjd(year ,month, day, hour, min, sec);
                MJD = JD1 - 2400000.5 + JD2;
                MJD_day = floor(MJD);
                sod = round((MJD - MJD_day) * 86400);
                satnum = str2double(str(34:35));
                obsdata.epoch(i, :) = [MJD_day sod];
                for sati = 1:satnum
                    str = fgets(fid);
                    % ------ GPS ------ %
                    if strcmp(str(1),'G') && OG_on
                        try
                            PRN = str2double(str(2:3));
                            obsdata.P1(i, PRN) = str2double(str(16*obsdata.OG(1)-12 : 16*obsdata.OG(1)+1));
                            obsdata.L1(i, PRN) = str2double(str(16*obsdata.OG(2)-12 : 16*obsdata.OG(2)+1));
                            obsdata.P2(i, PRN) = str2double(str(16*obsdata.OG(3)-12 : 16*obsdata.OG(3)+1));
                            obsdata.L2(i, PRN) = str2double(str(16*obsdata.OG(4)-12 : 16*obsdata.OG(4)+1));
                        catch
                            continue;
                        end
                    
                    % ------ BDS ------ %
                    elseif strcmp(str(1),'C') && OC_on
                        PRN = str2double(str(2:3));
                        obsdata.P1(i, PRN) = str2double(str(16*obsdata.OC(1)-12 : 16*obsdata.OC(1)+1));
                        obsdata.L1(i, PRN) = str2double(str(16*obsdata.OC(2)-12 : 16*obsdata.OC(2)+1));
                        obsdata.P2(i, PRN) = str2double(str(16*obsdata.OC(3)-12 : 16*obsdata.OC(3)+1));
                        obsdata.L2(i, PRN) = str2double(str(16*obsdata.OC(4)-12 : 16*obsdata.OC(4)+1));
                    
                    % ------ GAL ------ %
                    elseif strcmp(str(1),'E') && OE_on
                        PRN = str2double(str(2:3));
                        obsdata.P1(i, PRN) = str2double(str(16*obsdata.OE(1)-12 : 16*obsdata.OE(1)+1));
                        obsdata.L1(i, PRN) = str2double(str(16*obsdata.OE(2)-12 : 16*obsdata.OE(2)+1));
                        obsdata.P2(i, PRN) = str2double(str(16*obsdata.OE(3)-12 : 16*obsdata.OE(3)+1));
                        obsdata.L2(i, PRN) = str2double(str(16*obsdata.OE(4)-12 : 16*obsdata.OE(4)+1));
                    end
                end
            end
        end
        % 寻找并去除全为NaN值的行
        rowIsAllNaN = all(isnan(obsdata.P1), 2) & all(isnan(obsdata.P2), 2) & all(isnan(obsdata.L1), 2) & all(isnan(obsdata.L2), 2);
        obsdata.epoch(rowIsAllNaN, :) = [];
        obsdata.P1(rowIsAllNaN, :) = [];
        obsdata.P2(rowIsAllNaN, :) = [];
        obsdata.L1(rowIsAllNaN, :) = [];
        obsdata.L2(rowIsAllNaN, :) = [];
    else
        obsdata = -1;
    end
end