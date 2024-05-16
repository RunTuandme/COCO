function [dpsipr, depspr] = iauPr00(date1, date2)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
   DAS2R = 4.848136811095359935899141e-6;
   DJ00 = 2451545.0;
   DJC = 36525.0;

   %Precession and obliquity corrections (radians per century)
   PRECOR = -0.29965 * DAS2R;
   OBLCOR = -0.02524 * DAS2R;

   % Interval between fundamental epoch J2000.0 and given date (JC).
   t = ((date1 - DJ00) + date2) / DJC;

   % Precession rate contributions with respect to IAU 1976/80.
   dpsipr = PRECOR * t;
   depspr = OBLCOR * t;
end