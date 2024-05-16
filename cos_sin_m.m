function [c, s] = cos_sin_m(lam, order)
% DESCRIPTION:     Calculate cos(lam) and sin(lam).
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           lam, geocentric longitude.
%                  order, order of the gravitational field model.    
    c = zeros(1,order);
    s = zeros(1,order);
    for m = 1:order
        if m == 1
            c(m) = cos(lam);
            s(m) = sin(lam);
            continue
        end
        c(m) = cos(lam) * c(m-1) - sin(lam) * s(m-1);
        s(m) = cos(lam) * s(m-1) + sin(lam) * c(m-1);
%         c(m) = cos(m*lam);
%         s(m) = sin(m*lam);
    end
end