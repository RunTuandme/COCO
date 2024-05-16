function sp3data = HANDLE_readsp3(filename)
% DESCRIPTION:     Read the SP3 file.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    global obscodes
    fid = fopen(filename, 'r');
    str = fgets(fid);
    N_eops = str2double(str(33:39));
    sp3data.epoch = NaN(N_eops, 2);
    sp3data.pos = NaN(N_eops, 3, obscodes.num);
    sp3data.available = zeros(1, obscodes.num);
    i = 0;
    while ~feof(fid)
        str = fgets(fid);
        if strcmp(str(1), '*')
            i = i + 1;
            year = str2double(str(4:7));
            month = str2double(str(9:10));
            day = str2double(str(12:13));
            hour = str2double(str(15:16));
            min = str2double(str(18:19));
            sec = str2double(str(21:31));
            [JD1, JD2] = TOOL_countjd(year ,month, day, hour, min, sec);
            MJD = JD1 - 2400000.5 + JD2;
            MJD_day = floor(MJD);
            sod = round((MJD - MJD_day) * 86400);
            sp3data.epoch(i,:) = [MJD_day sod];
        end
        if strcmp(str(1:2), 'PG')
            PRN = str2double(str(3:4));
            X = str2double(str(6:18)) * 1e3;
            Y = str2double(str(20:32)) * 1e3;
            Z = str2double(str(34:46)) * 1e3;
            sp3data.pos(i, :, PRN) = [X Y Z];
            sp3data.available(PRN) = PRN;
        end 
    end
    % 寻找并去除全为NaN值的行
    sp3data.epoch = TOOL_clearNaN(sp3data.epoch);
    rowsWithNaN = any(all(isnan(sp3data.pos), 3), 2);
    sp3data.pos = sp3data.pos(~rowsWithNaN, :, :);
end