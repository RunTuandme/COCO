function [elong, phi, height] = pv2ob(x, y, z)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0


    % Given: xyz in WGS84
    % Returned:
    %**     elong        longitude (radians, east + ve)
    %**     phi          latitude (geodetic, radians)
    %**     height       height above ellipsoid (geodetic)
    
    a = 6378137.0;              % WGS84 equatorial radius
    f = 1.0 / 298.257223563;    % WGS84 flattening
    
    % Functions of ellipsoid parameters (with further validation of f).
    aeps2 = a*a * 1e-32;
    e2 = (2.0 - f) * f;
    e4t = e2*e2 * 1.5;
    ec2 = 1.0 - e2;
    
    ec = sqrt(ec2);
    b = a * ec;
    
    % Distance from polar axis squared.
    p2 = x*x + y*y;

    % Longitude.
    if p2 > 0.0
        elong = atan2(y, x);
    else
        elong = 0.0;
    end

    % Unsigned z-coordinate. 
   absz = abs(z);

    % Proceed unless polar case. 
   if ( p2 > aeps2 ) 

    % Distance from polar axis. 
      p = sqrt(p2);

    % Normalization. 
      s0 = absz / a;
      pn = p / a;
      zc = ec * s0;

    % Prepare Newton correction factors. 
      c0 = ec * pn;
      c02 = c0 * c0;
      c03 = c02 * c0;
      s02 = s0 * s0;
      s03 = s02 * s0;
      a02 = c02 + s02;
      a0 = sqrt(a02);
      a03 = a02 * a0;
      d0 = zc*a03 + e2*s03;
      f0 = pn*a03 - e2*c03;

    % Prepare Halley correction factor. 
      b0 = e4t * s02 * c02 * pn * (a0 - ec);
      s1 = d0*f0 - b0*s0;
      cc = ec * (f0*f0 - b0*c0);

    % Evaluate latitude and height. 
      phi = atan(s1/cc);
      s12 = s1 * s1;
      cc2 = cc * cc;
      height = (p*cc + absz*s1 - a * sqrt(ec2*s12 + cc2)) / sqrt(s12 + cc2);
   else 

    % Exception: pole. 
      phi = pi / 2.0;
      height = absz - b;
   end

    % Restore sign of latitude. 
   if  z < 0  
       phi = -phi;
   end
   
   elong = elong / pi * 180;
   phi = phi / pi * 180;