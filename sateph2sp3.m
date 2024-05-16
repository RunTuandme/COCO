function sateph2sp3(jdstart, sateph, step)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           jdstart, Julian Date of first epoch
%                  sateph, line 1-6 position(m) velosity(m/s)
%                  step, step size (seconds)
    N = size(sateph, 1);
    jdend = jdstart + (N-1) * step / 86400;
    jdlist = (jdstart: step/86400: jdend);
    mjdstart = jdstart - 2400000.5;
    [head_year, head_month, head_day, head_hour, head_minute, head_second] = jd2ce(jdstart);
    [head_year, head_month, head_day, head_hour, head_minute, head_second] = ...
        invert_date(head_year, head_month, head_day, head_hour, head_minute, head_second);
    [GPSweek, GPSsec] = mjd2gpst(mjdstart);
    
    fid = fopen(['./output/','shao',num2str(GPSweek),num2str(floor(GPSsec/86400)),'.sp3'], 'wt');
    
    fprintf(fid, '#cV%4d %2d %2d %2d %2d %11.8f %7d ORBIT CTS   SPP SHAO\n', ...
        head_year, head_month, head_day, head_hour, head_minute, head_second, N);
    fprintf(fid, '## %4d %15.8f %14.8f %5d %15.13f\n', ...
        GPSweek, GPSsec, step, floor(mjdstart), mjdstart-floor(mjdstart));
    fprintf(fid, '+    1   L01  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n');
    fprintf(fid, '+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n');
    fprintf(fid, '+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n');
    fprintf(fid, '+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n');
    fprintf(fid, '+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n');
    fprintf(fid, '++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n');
    fprintf(fid, '++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n');
    fprintf(fid, '++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n');
    fprintf(fid, '++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n');
    fprintf(fid, '++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n');
    fprintf(fid, '%%c L  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n');
    fprintf(fid, '%%f  1.2500000  1.025000000  0.00000000000  0.000000000000000\n');
    fprintf(fid, '%%f  0.0000000  0.000000000  0.00000000000  0.000000000000000\n');
    fprintf(fid, '%%i    0    0    0    0      0      0      0      0         0\n');
    fprintf(fid, '%%i    0    0    0    0      0      0      0      0         0\n');
    fprintf(fid, '/*                                                          \n');
    fprintf(fid, '/*                                                          \n');
    fprintf(fid, '/*                                                          \n');
    fprintf(fid, '/*Shanghai Astronomical Observatory LEO Orbit               \n');
    for i = 1: N
        jd = jdlist(i);
        [year, month, day, hour, minut, sec]=jd2ce(jd);
        [year, month, day, hour, minut, sec] = invert_date(year, month, day, hour, minut, sec);
        x = sateph(i, 1) * 1e-3;
        y = sateph(i, 2) * 1e-3;
        z = sateph(i, 3) * 1e-3;
        vx = sateph(i, 4) * 10;
        vy = sateph(i, 5) * 10;
        vz = sateph(i, 6) * 10;
        fprintf(fid, '*  %4d %2d %2d %2d %2d %11.8f\n', ...
            year, month, day, hour, minut, sec);
        fprintf(fid, 'PL01%14.6f%14.6f%14.6f%14.6f\n', ...
            x, y, z, 0);
        fprintf(fid, 'VL01%14.6f%14.6f%14.6f%14.6f\n', ...
            vx, vy, vz, 0);
    end
    fprintf(fid, 'EOF');
end

function [Y, M, D, H, m, S] = invert_date(Y, M, D, H, m, S)
    % dt: datetime
    dt = datetime([Y M D H m S]);
    DT = datetime(round(datenum(dt).*86400)./86400,'ConvertFrom','datenum');
    Y = DT.Year;
    M = DT.Month;
    D = DT.Day;
    H = DT.Hour;
    m = DT.Minute;
    S = DT.Second;
end