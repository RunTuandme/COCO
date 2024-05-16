function [rb, rp, rbp] = iauBp00(date1, date2)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
   % J2000.0 obliquity (Lieske et al. 1977)
   DJ00 = 2451545.0;
   DJC = 36525.0;
   DAS2R = 4.848136811095359935899141e-6;
   EPS0 = 84381.448 * DAS2R;

   % Interval between fundamental epoch J2000.0 and current date (JC).
   t = ((date1 - DJ00) + date2) / DJC;

   % Frame bias.
   [dpsibi, depsbi, dra0] = iauBi00;

   % Precession angles (Lieske et al. 1977) 
   psia77 = (5038.7784 + (-1.07259 + (-0.001147) * t) * t) * t * DAS2R;
   oma77  =       EPS0 + ((0.05127 + (-0.007726) * t) * t) * t * DAS2R;
   chia   = (  10.5526 + (-2.38064 + (-0.001125) * t) * t) * t * DAS2R;

   % Apply IAU 2000 precession corrections.
   [dpsipr, depspr] = iauPr00(date1, date2);
   psia = psia77 + dpsipr;
   oma  = oma77  + depspr;

   % Frame bias matrix: GCRS to J2000.0.
   rb = iauRx(-depsbi) * iauRy(dpsibi*sin(EPS0)) * iauRz(dra0);

   % Precession matrix: J2000.0 to mean of date. 
   rp = iauRz(chia) * iauRx(-oma) * iauRz(-psia) * iauRx(EPS0);

   % Bias-precession matrix: GCRS to mean of date.
   rbp = rp * rb;
end