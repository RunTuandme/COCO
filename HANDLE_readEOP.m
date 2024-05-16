function EOPmat = HANDLE_readEOP(fn)
% DESCRIPTION:     Read the EOP file.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    EOPmat = NaN(2e2, 7);
    fid = fopen(fn, 'r');
    i = 1;
    while ~feof(fid)
        T = fgetl(fid);
        mjd = round(str2double(T(8:15)));
        xp = str2double(T(20:27));
        yp = str2double(T(39:46));
        UT1_UTC = str2double(T(59:68));
        LOD = str2double(T(80:86));
        dX = str2double(T(97:106));
        dY = str2double(T(117:125));
        EOPmat(i,:) = [mjd xp yp UT1_UTC LOD dX dY];
        i = i + 1;
    end
end