function [sp3data, clkdata, obsdata] = HANDLE_interp(sp3data, clkdata, obsdata)
% DESCRIPTION:     Interpolate the data.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    global interpSP3posF
    global interpSP3velF
    global interpCLKF
    global obs_stepwithfile
    
    aclk = clkdata.clk(:,clkdata.available(clkdata.available~=0));
    clk_line = any(isnan(aclk), 2);
    clkdata.clk = clkdata.clk(clk_line==0,:);
    clkdata.epoch = clkdata.epoch(clk_line==0,:);

    % 确定时间基准 UTC min(MJD) 00:00:00
    MIN_MJD = min(sp3data.epoch(1,1), min(clkdata.epoch(1,1), obsdata.epoch(1,1)));
    sp3Day = sp3data.epoch(:,1) - MIN_MJD;
    clkDay = clkdata.epoch(:,1) - MIN_MJD;
    obsDay = obsdata.epoch(:,1) - MIN_MJD;
    sp3tlist = sp3Day * 86400 + sp3data.epoch(:,2);
    clktlist = clkDay * 86400 + clkdata.epoch(:,2);
    obstlist = obsDay * 86400 + obsdata.epoch(:,2);
    sp3data.epoch(:,3) = sp3tlist;
    clkdata.epoch(:,3) = clktlist;
    obsdata.epoch(:,3) = obstlist;
    
    % 重采样
    newTimeSeries = 1 : obs_stepwithfile : size(obstlist, 1);
    obstlist = obstlist(newTimeSeries);
    obsdata.epoch = obsdata.epoch(newTimeSeries,:);
    obsdata.P1 = obsdata.P1(newTimeSeries,:);
    obsdata.P2 = obsdata.P2(newTimeSeries,:);
    obsdata.L1 = obsdata.L1(newTimeSeries,:);
    obsdata.L2 = obsdata.L2(newTimeSeries,:);  
    
    % 去重
    [~, unique_rows_idx] = unique(sp3tlist, 'rows');
    sp3tlist = sp3tlist(unique_rows_idx);
    sp3data.pos = sp3data.pos(unique_rows_idx, :, :);
    
     % 星历插值梯度场
    [row, col, num] = size(sp3data.pos);
    sp3data.vel = NaN(row, col, num);
    
    timeinter = (obstlist(1):obstlist(end));
    sp3inter = NaN(length(timeinter), 3, num);
    for numsat = 1:num
        xlist = interp_larange_list(sp3tlist, sp3data.pos(:,1,numsat), timeinter);
        ylist = interp_larange_list(sp3tlist, sp3data.pos(:,2,numsat), timeinter);
        zlist = interp_larange_list(sp3tlist, sp3data.pos(:,3,numsat), timeinter);
        satn = [xlist ylist zlist];
        sp3inter(:,:,numsat) = satn;
    end
    
    interpSP3posF = griddedInterpolant(timeinter, sp3inter, 'spline');
    sp3data.vel = (interpSP3posF(sp3tlist+0.01) - sp3data.pos) / 0.01;

%     interpSP3posF = griddedInterpolant(sp3tlist, sp3data.pos, 'spline');
%     sp3data.vel = (interpSP3posF(sp3tlist+0.01) - sp3data.pos) / 0.01;
    
    % 对齐 data
    overlapping_arcs = (obsdata.epoch(:,3) >= sp3tlist(1)) & (obsdata.epoch(:,3) <= sp3tlist(end));
    obsdata.epoch = obsdata.epoch(overlapping_arcs, :);
    obsdata.P1 = obsdata.P1(overlapping_arcs, :);
    obsdata.P2 = obsdata.P2(overlapping_arcs, :);
    obsdata.L1 = obsdata.L1(overlapping_arcs, :);
    obsdata.L2 = obsdata.L2(overlapping_arcs, :);
    obstlist = (obsdata.epoch(:,3):obsdata.epoch(end,3));
    
    colHasNaN = any(isnan(clkdata.clk), 1);
    layerHasNaN = any(isnan(sp3data.pos), [1, 2]);
    NaNlogical = ( colHasNaN(:) | layerHasNaN(:) );
    layerHasNaN(1,1,:) = NaNlogical;
    sp3data.pos(:,:,layerHasNaN(1,1,:)) = NaN;
    clkdata.clk(:,NaNlogical') = NaN;
    
    interpCLKF = griddedInterpolant(clktlist, clkdata.clk);
    clkdata.clk = interpCLKF(obstlist);

    vector_pos = interpSP3posF(obstlist);

    sp3data.pos = vector_pos;
    interpSP3velF = griddedInterpolant(sp3tlist, sp3data.vel, 'spline');
    vector_vel = interpSP3velF(obstlist);
    sp3data.vel = vector_vel;
    
    % 对齐 epoch
    sp3data.epoch = obsdata.epoch;
    clkdata.epoch = obsdata.epoch;
    
    % 去掉观测值文件中频点丢失的数据
    nanIndices = isnan(obsdata.L1) | isnan(obsdata.L2) | isnan(obsdata.P1) | isnan(obsdata.P2);
    obsdata.L1(nanIndices) = NaN;
    obsdata.L2(nanIndices) = NaN;
    obsdata.P1(nanIndices) = NaN;
    obsdata.P2(nanIndices) = NaN;
    
    % 计算采样间隔
    obsdata.num = size(obsdata.epoch, 1);
    steplist = diff(obsdata.epoch(:,3));
    obsdata.step = mode(steplist);
    obsdata.epochnum = zeros(obsdata.num, 1);
    obsdata.epochnum(1) = 1;
    for i = 2: obsdata.num
        obsdata.epochnum(i) = obsdata.epochnum(i-1) + steplist(i-1) / obsdata.step;
    end
    
    %
end