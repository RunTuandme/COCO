function gmst = iauGmst00(uta, utb, tta, ttb)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
   DJ00 = 2451545.0;
   DJC = 36525.0;
   DAS2R = 4.848136811095359935899141e-6;

   % TT Julian centuries since J2000.0.
   t = ((tta - DJ00) + ttb) / DJC;

   % Greenwich Mean Sidereal Time, IAU 2000. 
   gmst = iauAnp(iauEra00(uta, utb) + ...
                   (     0.014506   + ...
                   (  4612.15739966 + ...
                   (     1.39667721 + ... 
                   (    -0.00009344 + ...
                   (     0.00001882 ) ...
          * t) * t) * t) * t) * DAS2R);
end