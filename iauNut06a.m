function [dpsi, deps] = iauNut06a(date1, date2)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    DJ00 = 2451545.0;
    DJC = 36525.0;
    
    % Interval between fundamental date J2000.0 and given date (JC). 
    t = ((date1 - DJ00) + date2) / DJC;
    
    % Factor correcting for secular variation of J2. 
    fj2 = -2.7774e-6 * t;
    
    % Obtain IAU 2000A nutation. 
    [dp, de] = iauNut00a_mex(date1, date2);
    
    % Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5).
    dpsi = dp + dp * (0.4697e-6 + fj2);
    deps = de + de * fj2;
end