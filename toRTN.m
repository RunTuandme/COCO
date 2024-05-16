function vector_rtn = toRTN(rx, ry, rz, vx, vy, vz, vector_GCRF)
% DESCRIPTION:     Convert vector from GCRF to RTN.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           rx, ry, rz, vx, vy, vz, reference satellite position and velocity in GCRF.
%                  vector_GCRF, vector in GCRF.
% OUTPUT:          vector_rtn, vector in RTN.

vector_r = [rx; ry; rz];
vector_v = [vx; vy; vz];

% 法向量
vector_n = cross(vector_r, vector_v);

% RTN 基向量
er = vector_r / norm(vector_r);
en = vector_n / norm(vector_n);
et = cross(en, er);

% 投影矩阵
B = [er'; et'; en' ];

vector_rtn = B * vector_GCRF;
end