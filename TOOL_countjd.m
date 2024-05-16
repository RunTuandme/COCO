function [jd1, jd2] = TOOL_countjd(year, month, day, hour, minut, sec)
% DESCRIPTION:     Count the Julian Day.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    if month == 1 || month == 2
        year = year - 1;
        month = month + 12;
    end 
    B = 2 - floor(year/100) + floor((0.25 * floor(year/100)));
    C = (sec/60 + minut) / 60 + hour;
    C = C/24;
    jd1 = floor(365.25 * (year + 4716)) + floor(30.6 * (month + 1)) + day + B - 1524.5;
    jd2 = C;
end

