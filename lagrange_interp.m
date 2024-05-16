
function Vq = lagrange_interp(X, V, Xq, n)
% DESCRIPTION:     Lagrange interpolation.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    N = length(X);
    if n > N
        error('阶数超过了数据点的数量');
    end

    Vq = 0;
    for i = 1:N
        L = 1;
        for j = 1:N
            if j ~= i
                L = L * (Xq - X(j)) / (X(i) - X(j));
            end
        end
        if i <= n
            Vq = Vq + V(i) * L;
        end
    end
end