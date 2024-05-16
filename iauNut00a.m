function [ddp, dde] = iauNut00a(date1, date2)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    DJ00 = 2451545.0;
    DJC = 36525.0;
    % Arcseconds to radians
    DAS2R = 4.848136811095359935899141e-6;
    % Arcseconds in a full circle
    TURNAS = 1296000.0;
    D2PI = 2 * pi;
    
    % Units of 0.1 microarcsecond to radians
    U2R = DAS2R / 1e7;
    
    % Luni-Solar nutation model
    % load('xls.mat', 'xls');
    global xls
    % Number of terms in the luni-solar nutation model
    NLS = 678;
    
    % Planetary nutation model
    % load('xpl.mat', 'xpl');
    global xpl
    NPL = 687;
    
    % Interval between fundamental date J2000.0 and given date (JC). 
    t = ((date1 - DJ00) + date2) / DJC;
    
    %% LUNI-SOLAR NUTATION
    
    % Mean anomaly of the Moon (IERS 2003).
    el = rem(485868.249036  + t * ( 1717915923.2178 + t * (31.8792 + t * (0.051635 + ...
             t * (- 0.00024470 ) ) ) ), TURNAS ) * DAS2R;
         
    % Mean anomaly of the Sun (MHB2000).
    elp = rem(1287104.79305  + ...
            t * (129596581.0481  + ...
            t * (-0.5532  + ...
            t * (0.000136  + ...
            t * (-0.00001149)))), TURNAS) * DAS2R;
        
    % Mean longitude of the Moon minus that of the ascending node (IERS 2003).
    f = rem(335779.526232 + ...
             t * (1739527262.8478 + ...
             t * (-12.7512 + ... 
             t * (-0.001037 + ...
             t * (0.00000417 ) ) ) ), TURNAS ) * DAS2R;
      
    % Mean elongation of the Moon from the Sun (MHB2000).
    d = rem(1072260.70369  + ...
          t * (1602961601.2090  + ... 
          t * (-6.3706  + ...
          t * (0.006593  + ...
          t * (-0.00003169)))), TURNAS) * DAS2R;
      
    % Mean longitude of the ascending node of the Moon (IERS 2003).
    om = rem(450160.398036 + ...
             t * (-6962890.5431 + ...
             t * (7.4722 + ...
             t * (0.007702 + ...
             t * (- 0.00005939 ) ) ) ), TURNAS ) * DAS2R;
         
    % Initialize the nutation values.
    dp = 0.0;
    de = 0.0;
    
    % Summation of luni-solar nutation series (in reverse order). 
   for i = 1:NLS

        % Argument and functions.
      arg =  rem(xls(i,1)  * el + ...
                 xls(i,2)  * elp + ...
                 xls(i,3)  * f + ...
                 xls(i,4)  * d + ...
                 xls(i,5)  * om, D2PI);
      sarg = sin(arg);
      carg = cos(arg);

      % Term. 
      dp = dp + (xls(i,6) + xls(i,7) * t) * sarg + xls(i,8) * carg;
      de = de + (xls(i,9) + xls(i,10) * t) * carg + xls(i,11) * sarg;
   end
   
   % Convert from 0.1 microarcsec units to radians.
   dpsils = dp * U2R;
   depsls = de * U2R;
   
   %% PLANETARY NUTATION
   % Mean anomaly of the Moon (MHB2000).
   al = rem(2.35555598 + 8328.6914269554 * t, D2PI);
   
   % Mean longitude of the Moon minus that of the ascending node (MHB2000).
   af = rem(1.627905234 + 8433.466158131 * t, D2PI);
   
   % Mean elongation of the Moon from the Sun (MHB2000).
   ad = rem(5.198466741 + 7771.3771468121 * t, D2PI);
   
   % Mean longitude of the ascending node of the Moon (MHB2000).
   aom = rem(2.18243920 - 33.757045 * t, D2PI);
   
   % General accumulated precession in longitude (IERS 2003).
   apa = (0.024381750 + 0.00000538691 * t) * t;
   
   % Planetary longitudes, Mercury through Uranus (IERS 2003).
   alme = rem(4.402608842 + 2608.7903141574 * t, D2PI);
   alve = rem(3.176146697 + 1021.3285546211 * t, D2PI);
   alea = rem(1.753470314 + 628.3075849991 * t, D2PI);
   alma = rem(6.203480913 + 334.0612426700 * t, D2PI);
   alju = rem(0.599546497 + 52.9690962641 * t, D2PI);
   alsa = rem(0.874016757 + 21.3299104960 * t, D2PI);
   alur = rem(5.481293872 + 7.4781598567 * t, D2PI);
   
   % Neptune longitude (MHB2000).
   alne = rem(5.321159000 + 3.8127774000 * t, D2PI);
   
   % Initialize the nutation values.
   dp = 0.0;
   de = 0.0;
   
   % Summation of planetary nutation series (in reverse order).
   for i = 1: NPL
      % Argument and functions.
      arg =  rem(xpl(i,1) * al   + ...
                 xpl(i,2) * af   + ...
                 xpl(i,3) * ad   + ...
                 xpl(i,4) * aom  + ...
                 xpl(i,5) * alme + ...
                 xpl(i,6) * alve + ...
                 xpl(i,7) * alea + ...
                 xpl(i,8) * alma + ...
                 xpl(i,9) * alju + ...
                 xpl(i,10) * alsa + ...
                 xpl(i,11) * alur + ...
                 xpl(i,12) * alne + ...
                 xpl(i,13) * apa, D2PI);
      sarg = sin(arg);
      carg = cos(arg);

      % Term
      dp = dp + xpl(i,14) * sarg + xpl(i,15) * carg;
      de = de + xpl(i,16) * sarg + xpl(i,17) * carg;
   end
   
   % Convert from 0.1 microarcsec units to radians.
   dpsipl = dp * U2R;
   depspl = de * U2R;
   
   % Add luni-solar and planetary components.
   ddp = dpsils + dpsipl;
   dde = depsls + depspl;
end

