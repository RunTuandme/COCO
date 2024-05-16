function s = iauS06(date1, date2, x, y)
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
    
    % Interval between fundamental date J2000.0 and given date (JC). 
    t = ((date1 - DJ00) + date2) / DJC;

    % Terms of order t^0
    % load('TERM.mat', 'TERM');
    global TERM
    
    % Terms of order t^1
    TERM1 = ...
    [ 0,  0,  0,  0,  2,  0,  0,  0,    -0.07e-6,   3.57e-6 
      0,  0,  0,  0,  1,  0,  0,  0,     1.73e-6,  -0.03e-6 
      0,  0,  2, -2,  3,  0,  0,  0,     0.00e-6,   0.48e-6 ];
  
    % Terms of order t^2
    TERM2 = [
      0,  0,  0,  0,  1,  0,  0,  0,   743.52e-6,  -0.17e-6
      0,  0,  2, -2,  2,  0,  0,  0,    56.91e-6,   0.06e-6
      0,  0,  2,  0,  2,  0,  0,  0,     9.84e-6,  -0.01e-6
      0,  0,  0,  0,  2,  0,  0,  0,    -8.85e-6,   0.01e-6
      0,  1,  0,  0,  0,  0,  0,  0,    -6.38e-6,  -0.05e-6
      1,  0,  0,  0,  0,  0,  0,  0,    -3.07e-6,   0.00e-6
      0,  1,  2, -2,  2,  0,  0,  0,     2.23e-6,   0.00e-6
      0,  0,  2,  0,  1,  0,  0,  0,     1.67e-6,   0.00e-6
      1,  0,  2,  0,  2,  0,  0,  0,     1.30e-6,   0.00e-6
      0,  1, -2,  2, -2,  0,  0,  0,     0.93e-6,   0.00e-6
      
      
      1,  0,  0, -2,  0,  0,  0,  0,     0.68e-6,   0.00e-6
      0,  0,  2, -2,  1,  0,  0,  0,    -0.55e-6,   0.00e-6
      1,  0, -2,  0, -2,  0,  0,  0,     0.53e-6,   0.00e-6
      0,  0,  0,  2,  0,  0,  0,  0,    -0.27e-6,   0.00e-6
      1,  0,  0,  0,  1,  0,  0,  0,    -0.27e-6,   0.00e-6
      1,  0, -2, -2, -2,  0,  0,  0,    -0.26e-6,   0.00e-6
      1,  0,  0,  0, -1,  0,  0,  0,    -0.25e-6,   0.00e-6
      1,  0,  2,  0,  1,  0,  0,  0,     0.22e-6,   0.00e-6
      2,  0,  0, -2,  0,  0,  0,  0,    -0.21e-6,   0.00e-6
      2,  0, -2,  0, -1,  0,  0,  0,     0.20e-6,   0.00e-6
      
      
      0,  0,  2,  2,  2,  0,  0,  0,     0.17e-6,   0.00e-6
      2,  0,  2,  0,  2,  0,  0,  0,     0.13e-6,   0.00e-6
      2,  0,  0,  0,  0,  0,  0,  0,    -0.13e-6,   0.00e-6
      1,  0,  2, -2,  2,  0,  0,  0,    -0.12e-6,   0.00e-6
      0,  0,  2,  0,  0,  0,  0,  0,    -0.11e-6,   0.00e-6 ];

    TERM3 = ...
    [ 0,  0,  0,  0,  1,  0,  0,  0,     0.30e-6, -23.42e-6 
      0,  0,  2, -2,  2,  0,  0,  0,    -0.03e-6,  -1.46e-6 
      0,  0,  2,  0,  2,  0,  0,  0,    -0.01e-6,  -0.25e-6 
      0,  0,  0,  0,  2,  0,  0,  0,     0.00e-6,   0.23e-6 ];
  
    TERM4 = ...
    [ 0,  0,  0,  0,  1,  0,  0,  0,    -0.26e-6,  -0.01e-6 ];
    
    % Polynomial coefficients
    sp = [94.00e-6,3808.65e-6,-122.68e-6,-72574.11e-6,27.98e-6,15.62e-6];
    
    % Mean anomaly of the Moon.
    fa(1) = rem(           485868.249036  + ...
             t * ( 1717915923.2178 + ...
             t * (         31.8792 + ...
             t * (          0.051635 + ...
             t * (        - 0.00024470 ) ) ) ), TURNAS ) * DAS2R;

    % Mean anomaly of the Sun.
    fa(2) = rem(         1287104.793048 + ...
             t * ( 129596581.0481 + ...
             t * (       - 0.5532 + ...
             t * (         0.000136 + ...
             t * (       - 0.00001149 ) ) ) ), TURNAS ) * DAS2R;

    % Mean longitude of the Moon minus that of the ascending node.
    fa(3) = rem(           335779.526232 + ...
             t * ( 1739527262.8478 + ...
             t * (       - 12.7512 + ...
             t * (        - 0.001037 + ...
             t * (          0.00000417 ) ) ) ), TURNAS ) * DAS2R;

    % Mean elongation of the Moon from the Sun. 
    fa(4) = rem(          1072260.703692 + ...
             t * ( 1602961601.2090 + ...
             t * (        - 6.3706 + ...
             t * (          0.006593 + ...
             t * (        - 0.00003169 ) ) ) ), TURNAS ) * DAS2R;

    % Mean longitude of the ascending node of the Moon. 
    fa(5) = rem(          450160.398036 + ...
             t * ( - 6962890.5431 + ...
             t * (         7.4722 + ...
             t * (         0.007702 + ...
             t * (       - 0.00005939 ) ) ) ), TURNAS ) * DAS2R;

    % Mean longitude of Venus. 
    fa(6) = rem(3.176146697 + 1021.3285546211 * t, D2PI);

    % Mean longitude of Earth. 
    fa(7) = rem(1.753470314 + 628.3075849991 * t, D2PI);

    % General precession in longitude. 
    fa(8) = (0.024381750 + 0.00000538691 * t) * t;
   
    w0 = sp(1);
    w1 = sp(2);
    w2 = sp(3);
    w3 = sp(4);
    w4 = sp(5);
    w5 = sp(6);
   
    for i = 1:33
        a = 0.0;
        for j = 1:8
            a = a + TERM(i,j) * fa(j);
        end
        w0 = w0 + TERM(i,9) * sin(a) + TERM(i,10) * cos(a);
    end

    for i = 1:3
        a = 0.0;
        for j = 1:8
            a = a + TERM1(i,j) * fa(j);
        end
        w1 = w1 + TERM1(i,9) * sin(a) + TERM1(i,10) * cos(a);
    end
    
    for i = 1:25
        a = 0.0;
        for j = 1:8
            a = a + TERM2(i,j) * fa(j);
        end
        w2 = w2 + TERM2(i,9) * sin(a) + TERM2(i,10) * cos(a);
    end
    
    for i = 1:4
        a = 0.0;
        for j = 1:8
            a = a + TERM3(i,j) * fa(j);
        end
        w3 = w3 + TERM3(i,9) * sin(a) + TERM3(i,10) * cos(a);
    end

    for i = 1:1
        a = 0.0;
        for j = 1:8
            a = a + TERM4(i,j) * fa(j);
        end
        w4 = w4 + TERM4(i,9) * sin(a) + TERM4(i,10) * cos(a);
    end

   s = (w0 + ...
       (w1 + ...
       (w2 + ...
       (w3 + ...
       (w4 + ...
        w5 * t) * t) * t) * t) * t) * DAS2R - x*y/2.0;
end
