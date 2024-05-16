function eps0 = iauObl80(date1, date2)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
   DJ00 = 2451545.0;
   DJC = 36525.0;
   DAS2R = 4.848136811095359935899141e-6;
   
   % Interval between fundamental epoch J2000.0 and given date (JC).
   t = ((date1 - DJ00) + date2) / DJC;

   % Mean obliquity of date.
   eps0 = DAS2R * (84381.448  + ...
                  (-46.8150   + ...
                  (-0.00059   + ...
                  ( 0.001813) * t) * t) * t);             
end