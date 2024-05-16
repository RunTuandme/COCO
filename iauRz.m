function C = iauRz(theta)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
C = [cos(theta)     sin(theta)          0;
     -sin(theta)    cos(theta)          0;
        0               0               1];
end