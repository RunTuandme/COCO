function eect = iauEect00(date1, date2)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
   DJ00 = 2451545.0;
   DJC = 36525.0;
   TURNAS = 1296000.0;
   DAS2R = 4.848136811095359935899141e-6;
   D2PI = 2 * pi;
    
   % Terms of order t^0
   e0 = [0,  0,  0,  0,  1,  0,  0,  0, 2640.96e-6, -0.39e-6 ;
      0,  0,  0,  0,  2,  0,  0,  0,   63.52e-6, -0.02e-6 ;
      0,  0,  2, -2,  3,  0,  0,  0,   11.75e-6,  0.01e-6 ;
      0,  0,  2, -2,  1,  0,  0,  0,   11.21e-6,  0.01e-6 ;
      0,  0,  2, -2,  2,  0,  0,  0,   -4.55e-6,  0.00e-6 ;
      0,  0,  2,  0,  3,  0,  0,  0,    2.02e-6,  0.00e-6 ;
      0,  0,  2,  0,  1,  0,  0,  0,    1.98e-6,  0.00e-6 ;
      0,  0,  0,  0,  3,  0,  0,  0,   -1.72e-6,  0.00e-6 ;
      0,  1,  0,  0,  1,  0,  0,  0,   -1.41e-6, -0.01e-6 ;
      0,  1,  0,  0, -1,  0,  0,  0,   -1.26e-6, -0.01e-6 ;
      1,  0,  0,  0, -1,  0,  0,  0,   -0.63e-6,  0.00e-6 ;
      1,  0,  0,  0,  1,  0,  0,  0,   -0.63e-6,  0.00e-6 ;
      0,  1,  2, -2,  3,  0,  0,  0,    0.46e-6,  0.00e-6 ;
      0,  1,  2, -2,  1,  0,  0,  0,    0.45e-6,  0.00e-6 ;
      0,  0,  4, -4,  4,  0,  0,  0,    0.36e-6,  0.00e-6 ;
      0,  0,  1, -1,  1, -8, 12,  0,   -0.24e-6, -0.12e-6 ;
      0,  0,  2,  0,  0,  0,  0,  0,    0.32e-6,  0.00e-6 ;
      0,  0,  2,  0,  2,  0,  0,  0,    0.28e-6,  0.00e-6 ;
      1,  0,  2,  0,  3,  0,  0,  0,    0.27e-6,  0.00e-6 ;
      1,  0,  2,  0,  1,  0,  0,  0,    0.26e-6,  0.00e-6 ;
      0,  0,  2, -2,  0,  0,  0,  0,   -0.21e-6,  0.00e-6 ;
      0,  1, -2,  2, -3,  0,  0,  0,    0.19e-6,  0.00e-6 ;
      0,  1, -2,  2, -1,  0,  0,  0,    0.18e-6,  0.00e-6 ;
      0,  0,  0,  0,  0,  8,-13, -1,   -0.10e-6,  0.05e-6 ;
      0,  0,  0,  2,  0,  0,  0,  0,    0.15e-6,  0.00e-6 ;
      2,  0, -2,  0, -1,  0,  0,  0,   -0.14e-6,  0.00e-6 ;
      1,  0,  0, -2,  1,  0,  0,  0,    0.14e-6,  0.00e-6 ;
      0,  1,  2, -2,  2,  0,  0,  0,   -0.14e-6,  0.00e-6 ;
      1,  0,  0, -2, -1,  0,  0,  0,    0.14e-6,  0.00e-6 ;
      0,  0,  4, -2,  4,  0,  0,  0,    0.13e-6,  0.00e-6 ;
      0,  0,  2, -2,  4,  0,  0,  0,   -0.11e-6,  0.00e-6 ;
      1,  0, -2,  0, -3,  0,  0,  0,    0.11e-6,  0.00e-6 ;
      1,  0, -2,  0, -1,  0,  0,  0,    0.11e-6,  0.00e-6 ];

   % Terms of order t^1 */
   e1 = [ 0,  0,  0,  0,  1,  0,  0,  0,  -0.87e-6,  0.00e-6 ];

   % Number of terms in the series 
   NE0 = size(e0, 1);
   NE1 = size(e1, 1);

   % Interval between fundamental epoch J2000.0 and current date (JC).
   t = ((date1 - DJ00) + date2) / DJC;

   % Fundamental Arguments (from IERS Conventions 2003)

   % Mean anomaly of the Moon. 
   fa(1) = rem(485868.249036  + t * ( 1717915923.2178 + t * (31.8792 + t * (0.051635 + ...
             t * (- 0.00024470 ) ) ) ), TURNAS ) * DAS2R;

   % Mean anomaly of the Sun.
   fa(2) = rem(1287104.793048 + t * ( 129596581.0481 + t * (- 0.5532 + t * (0.000136 + ...
             t * (- 0.00001149 ) ) ) ), TURNAS ) * DAS2R;

   % Mean longitude of the Moon minus that of the ascending node. 
   fa(3) = rem(335779.526232 + t * ( 1739527262.8478 + t * (- 12.7512 + t * (- 0.001037 + ...
             t * ( 0.00000417 ) ) ) ), TURNAS ) * DAS2R;

   % Mean elongation of the Moon from the Sun.
   fa(4) = rem(1072260.703692 + t * ( 1602961601.2090 + t * (- 6.3706 + t * ( 0.006593 + ...
             t * ( - 0.00003169 ) ) ) ), TURNAS ) * DAS2R;

   % Mean longitude of the ascending node of the Moon.
   fa(5) = rem( 450160.398036 + t * ( - 6962890.5431 + t * ( 7.4722 + t * ( 0.007702 + ...
             t * ( - 0.00005939 ) ) ) ), TURNAS ) * DAS2R;

   % Mean longitude of Venus.
   fa(6) = rem(3.176146697 + 1021.3285546211 * t, D2PI);

   % Mean longitude of Earth. 
   fa(7) = rem(1.753470314 + 628.3075849991 * t, D2PI);

   % General precession in longitude.
   fa(8) = (0.024381750 + 0.00000538691 * t) * t;

   % Evaluate the EE complementary terms.
   s0 = 0.0;
   s1 = 0.0;

   for i = NE0: -1: 1
      a = 0.0;
      for j = 1:8
         a = a + e0(i,j) * fa(j);
      end
      s0 = s0 + e0(i,9) * sin(a) + e0(i,10) * cos(a);
   end

   for i = NE1: -1: 1
      a = 0.0;
      for j = 1:8
         a = a + e1(i,j) * fa(j);
      end
      s1 = s1 + e1(i,9) * sin(a) + e1(i,10) * cos(a);
   end

   eect = (s0 + s1 * t ) * DAS2R;
end