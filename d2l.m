function [lon, lat, rho] = d2l(x,y,z)
% DESCRIPTION:     x y z to geocentric longitude and latitude.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           x, y, z, ITRF position.
    rho = norm([x y z]);
    lat = abs(asin( z / rho ));
    if z < 0
        lat = -lat;
    end
    s = [x y 0];
    e = [1 0 0];
    lon = acos(dot(s,e) / norm(s));
    if y < 0
        lon = 2 * pi - lon;
    end
end