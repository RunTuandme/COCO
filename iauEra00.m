function theta = iauEra00(dj1, dj2)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    DJ00 = 2451545.0;
    D2PI = 2 * pi;
    
    % Days since fundamental epoch.
    if (dj1 < dj2) 
      d1 = dj1;
      d2 = dj2;
    else 
      d1 = dj2;
      d2 = dj1;
    end
    t = d1 + (d2- DJ00);
    
    % Fractional part of T (days).
    f = rem(d1, 1.0) + rem(d2, 1.0);
    
    % Earth rotation angle at this UT1.
    theta = iauAnp(D2PI * (f + 0.7790572732640 ...
                            + 0.00273781191135448 * t));          
end