function w = iauAnp(a)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    D2PI = 6.283185307179586476925287;
    w = rem(a, D2PI);
   if w < 0
       w = w + D2PI;
   end
end