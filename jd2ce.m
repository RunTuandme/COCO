function [year, month, day, hour, minut, sec] = jd2ce(jd)
% DESCRIPTION:     Julian Date to Calendar Date.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           jd, Julian Date.
    JDN = jd + 0.5;
    Z = floor(JDN);
    F = JDN - Z;
    if Z < 2299161
        A = Z;
    else
        a = floor((Z - 2305447.5) / 36524.25);
        A = Z + 10 + a - floor(a/4);
    end
    B = A + 1524;
    C = floor((B - 122.1) / 365.25);
    D = floor(365.25*C);
    E = floor((B - D) / 30.6001);
    day_m = B - D - floor(30.6001*E) + F;
    day = floor(day_m);
    if E < 14
        month = E - 1;
    elseif E < 16
        month = E - 13;
    end
    if month > 2
        year = C - 4716;
    elseif month == 1 || month == 2
        year = C - 4715;
    end
    G = (day_m - day) * 24;
    hour = floor(G);
    H = (G - hour) * 60;
    minut = floor(H);
    I = (H - minut) * 60;
    sec = I;
    
end