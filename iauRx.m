function C = iauRx(theta)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
C = [   1              0               0;
        0          cos(theta)      sin(theta);
        0          -sin(theta)     cos(theta)];
end