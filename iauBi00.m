function [dpsibi, depsbi, dra] = iauBi00
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
   % The frame bias corrections in longitude and obliquity
   DAS2R = 4.848136811095359935899141e-6;
   DPBIAS = -0.041775  * DAS2R;
   DEBIAS = -0.0068192 * DAS2R;

   % The ICRS RA of the J2000.0 equinox (Chapront et al., 2002)
   DRA0 = -0.0146 * DAS2R;

   % Return the results (which are fixed).
   dpsibi = DPBIAS;
   depsbi = DEBIAS;
   dra = DRA0;
end