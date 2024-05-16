function acc = experience(x, y, z, vx, vy, vz, CC, Cn, Sn)
% DESCRIPTION:     Calculate the acceleration of empirical force.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           x, y, z, vx, vy, vz, satellite position and velocity.
%                  CC, Cn, Sn, coefficients of the empirical force model.
    u = 3.9860044150E+14;
    vector_r = [x; y; z];
    vector_v = [vx; vy; vz];
    r = norm(vector_r);
    v = norm(vector_v);
    a = 2 / r - v^2 / u;
    a = 1 / a;
    ecosE = 1 - r / a;
    esinE =(x*vx + y*vy + z*vz)/sqrt(u*a);
    e = sqrt(ecosE^2 + esinE^2);
    cosf = (a-r) / (e*r) - (a*e) / r;
    sinf = a * sqrt(1-e^2) * esinE / (e * r); 
    % 法向量
    vector_n = cross(vector_r, vector_v);
    % RTN 基向量
    er = vector_r / r;
    et = vector_v / v;
    en = vector_n / norm(vector_n);
    % 投影矩阵(J2000-RTN)
    B = [er'; et'; en'];
    acc = B' * (CC + Cn * cosf + Sn * sinf) * 1e-9;
end