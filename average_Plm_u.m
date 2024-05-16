function [p, p_] = average_Plm_u(u, order)
% DESCRIPTION:     Calculate p and p_ of the Legendre polynomials.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           u, sine function of geocentric latitude.
%                  order, order of the gravitational field model.
    p = zeros(order+1,order+1);
    p_ = zeros(order,order);
    
    p(1,1) = 1.732050807568877 * sqrt(1-u^2);
    for l = 2:order+1
        p(l,l) = sqrt( (2*l+1)/(2*l) ) * sqrt(1-u^2) * p(l-1,l-1);
    end
    for l=2:order+1
        p(l,l-1) = sqrt((2*l+1)) * u * p(l-1,l-1);
    end
    for l = 3:order+1
        for m = 1:l-2
            p(l,m) = sqrt((2*l+1)*(2*l-1)/(l+m)/(l-m)) * u * p(l-1,m) -...
                sqrt((2*l+1)*(l-1+m)*(l-1-m)/(2*l-3)/(l+m)/(l-m)) * p(l-2,m);
        end
    end
    
    for l = 1:order
        for m = 1:l
            p_(l,m) = 1 / sqrt(1-u^2)...
                * (sqrt( (l+m+1)*(l-m) )...
                * p(l,m+1) - m * u / sqrt(1-u^2) * p(l,m));
        end
    end
end