function clk = TOOL_getCLK(t0, activeS)
% DESCRIPTION:     Get GNSS clock.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    % t0: 积秒
    global interpCLKF
    clk = interpCLKF(t0);
    clk = clk(:, activeS);
    clk = diag(clk);
end