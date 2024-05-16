function clkdata = HANDLE_readclk(filename)
% DESCRIPTION:     Read the CLK file.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    global obscodes
    fid = fopen(filename, 'r');
    N_eops = 92160;
    clkdata.epoch = NaN(N_eops, 2);
    clkdata.clk = NaN(N_eops, obscodes.num);
    clkdata.available = zeros(1, obscodes.num);
    i = 0;
    start_read = false;
    mjd_day_tmp = 0;
    sod_tmp = 0;
    COD_bias = 0;
    while ~feof(fid)
        str = fgets(fid);
        if contains(str, 'ANALYSIS CENTER')
            if contains(str, 'CODE')
                COD_bias = 5;
            end
            continue
        end
        if contains(str, 'END OF HEADER')
            start_read = true;
            continue;
        end
        if start_read
            if strcmp(str(1:4), 'AS G')
                PRN = str2double(str(5:6));
                year = str2double(str(9+COD_bias:12+COD_bias));
                month = str2double(str(14+COD_bias:15+COD_bias));
                day = str2double(str(17+COD_bias:18+COD_bias));
                hour = str2double(str(20+COD_bias:21+COD_bias));
                min = str2double(str(23+COD_bias:24+COD_bias));
                sec = str2double(str(26+COD_bias:34+COD_bias));
                clk = str2double(str(41+COD_bias:59+COD_bias));
                [JD1, JD2] = TOOL_countjd(year ,month, day, hour, min, sec);
                MJD = JD1 - 2400000.5 + JD2;
                MJD_day = floor(MJD);
                sod = round((MJD - MJD_day) * 86400);
                if ~(MJD_day == mjd_day_tmp && sod == sod_tmp)
                    mjd_day_tmp = MJD_day;
                    sod_tmp = sod;
                    i = i + 1;
                    clkdata.epoch(i, :) = [MJD_day sod];
                end    
                clkdata.clk(i, PRN) = clk;
                clkdata.available(PRN) = PRN;
            end  
        end
    end
    rowIsAllNaN = all(isnan(clkdata.clk), 2);
    clkdata.epoch(rowIsAllNaN, :) = [];
    clkdata.clk(rowIsAllNaN, :) = [];
end

% function clkdata = HANDLE_readclk(filename)
%     fid = fopen(filename, 'r');
%     N_eops = 92160;
%     clkdata.epoch = zeros(N_eops, 2);
%     clkdata.clk = zeros(N_eops, 105);
%     clkdata.available = zeros(1, 105);
%     i = 0;
%     start_read = false;
%     mjd_day_tmp = 0;
%     sod_tmp = 0;
%     COD_bias = 0;
%     
%     % Regular expression patterns for extracting values
%     asgPattern = '^AS G(\d{2})\s+(\d{4})(\d{2})(\d{2})\s+(\d{2})(\d{2})(\d{2}\.\d+)\s+([\d\s.]+)$';
%     codPattern = '^COD';
%     
%     while ~feof(fid)
%         str = fgets(fid);
%         if regexp(str, codPattern)
%             COD_bias = 5;
%             continue
%         end
%         if contains(str, 'END OF HEADER')
%             start_read = true;
%             continue;
%         end
%         if start_read
%             tokens = regexp(str, asgPattern, 'tokens');
%             if ~isempty(tokens)
%                 tokens = tokens{1};
%                 PRN = str2double(tokens{1});
%                 year = str2double(tokens{2});
%                 month = str2double(tokens{3});
%                 day = str2double(tokens{4});
%                 hour = str2double(tokens{5});
%                 min = str2double(tokens{6});
%                 sec = str2double(tokens{7});
%                 clk = str2double(tokens{8});
%                 jd = TOOL_countjd(year ,month, day, hour, min, sec, 0);
%                 mjd_day = floor(jd - 2400000.5);
%                 sod = round((jd - 2400000.5 - mjd_day) * 86400);
%                 if ~(mjd_day == mjd_day_tmp && sod == sod_tmp)
%                     mjd_day_tmp = mjd_day;
%                     sod_tmp = sod;
%                     i = i + 1;
%                     clkdata.epoch(i, :) = [mjd_day sod];
%                 end
%                 clkdata.clk(i, PRN) = clk;
%                 clkdata.available(PRN) = PRN;
%             end
%         end
%     end
%     
%     fclose(fid);
%     
%     % Remove rows with all NaN values
%     rowIsAllNaN = all(isnan(clkdata.clk), 2);
%     clkdata.epoch(rowIsAllNaN, :) = [];
%     clkdata.clk(rowIsAllNaN, :) = [];
% end