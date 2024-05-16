function [rel] = relativity(x, y, z, vx, vy, vz)
% DESCRIPTION:     Compute the relativity.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    beta = 1;
    gamma = 1;
    GM = 3.9860044150E+14;
    c = 2.99792458E+8;
    vector_r = [x; y; z];
    vector_v = [vx; vy; vz];
    r = norm(vector_r);
    v = norm(vector_v);
    
    A = (2 * (beta + gamma) * GM / r - gamma * v^2) * vector_v;
    B = 2 * (1 + gamma) * (vector_r' * vector_v) * vector_v;
    RL1 = GM / c^2 / r^3 * (A + B); % Schwarzschild 项
    RL2 = 0; % 测地岁差项
    RL3 = 0; % Lense - Thirring 岁差项
    rel = RL1 + RL2 + RL3;
end

