function [epsa, rb, rp, rbp, rn, rbpn] = iauPn00(date1, date2, dpsi, deps)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
   % IAU 2000 precession-rate adjustments.
   [~, depspr] = iauPr00(date1, date2);

   % Mean obliquity, consistent with IAU 2000 precession-nutation.
   epsa = iauObl80(date1, date2) + depspr;

   % Frame bias and precession matrices and their product.
   [rb, rp, rbp] = iauBp00(date1, date2);

   % Nutation matrix. 
   rn = iauRx(-(epsa + deps)) * iauRz(-dpsi) * iauRx(epsa);

   % Bias-precession-nutation matrix (classical).
   rbpn = rn * rbp;
end