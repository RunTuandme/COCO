function ac = thirdbody(x, y, z, JPL)
% DESCRIPTION:     Calculate the third body acceleration.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           x, y, z, satellite position.
%                  JPL, JPL results, the format is the gather of the output of the function JPL_Eph_DE440.m. e.g. [r_Mercury,r_Venus].
    r_sat = [x; y; z];
    r_Mercury = JPL(:,1);
    r_Venus = JPL(:,2);
%     r_Earth = JPL(:,3);
    r_Mars = JPL(:,4);
    r_Jupiter = JPL(:,5);
    r_Saturn = JPL(:,6);
    r_Uranus = JPL(:,7);
    r_Neptune = JPL(:,8);
    r_Pluto = JPL(:,9);
    r_Moon = JPL(:,10);
    r_Sun = JPL(:,11);
%     r_SunSSB = JPL(:,12);
%     gm_sun = 1.327124400179870e+20;
%     gm_moon = 4.902800582147763e+12;
    % Gravitational coefficients
    GM_Earth   = 398600.4415e9;         			   % [m^3/s^2]; GGM03C & GGM03S
    GM_Sun     = 132712440041.279419e9; 			   % [m^3/s^2]; DE440
    GM_Moon    = GM_Earth/81.3005682214972154;         % [m^3/s^2]; DE440
    GM_Mercury = 22031.868551e9; 		  			   % [m^3/s^2]; DE440
    GM_Venus   = 324858.592000e9;       			   % [m^3/s^2]; DE440
    GM_Mars    = 42828.375816e9;	      			   % [m^3/s^2]; DE440
    GM_Jupiter = 126712764.100000e9;    			   % [m^3/s^2]; DE440
    GM_Saturn  = 37940584.841800e9;     			   % [m^3/s^2]; DE440
    GM_Uranus  = 5794556.400000e9;      			   % [m^3/s^2]; DE440
    GM_Neptune = 6836527.100580e9;      			   % [m^3/s^2]; DE440
    GM_Pluto   = 975.500000e9;	      			       % [m^3/s^2]; DE440

    ac = N3tosat(r_sat, r_Moon, GM_Moon) + ...
         N3tosat(r_sat, r_Sun, GM_Sun) + ...
         N3tosat(r_sat, r_Mercury, GM_Mercury) + ...
         N3tosat(r_sat, r_Venus, GM_Venus) + ...
         N3tosat(r_sat, r_Mars, GM_Mars) + ...
         N3tosat(r_sat, r_Jupiter, GM_Jupiter) + ...
         N3tosat(r_sat, r_Saturn, GM_Saturn) + ...
         N3tosat(r_sat, r_Uranus, GM_Uranus) + ...
         N3tosat(r_sat, r_Neptune, GM_Neptune) + ...
         N3tosat(r_sat, r_Pluto, GM_Pluto);
%     ac = N3tosat(r_sat, r_Moon, GM_Moon) + ...
%           N3tosat(r_sat, r_Sun, GM_Sun);
end

function a = N3tosat(r, s, GM)
    d = r - s;
    a = -GM * ( d/(norm(d)^3) + s/(norm(s)^3) );
end
