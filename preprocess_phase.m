function [phasedata, pseudodata] = preprocess_phase(obsdata, sp3data, clkdata, detect_pos, spp_clk, start_t)

%  输出：phasedata 格式
% | 载波相位 | 历元号 | 整周未知数号 | 卫星号 | gnsssat_x | gnsssat_y | gnsssat_z | 导航星钟差 |
% 接收机钟差 | NaN | NaN | 整秒时 | 儒略日（整天） | 儒略日（小数天） 
%  输出：pseudodata 格式
% | 伪距 | 历元号 | NaN | 卫星号 | gnsssat_x | gnsssat_y | gnsssat_z | 导航星钟差 |
% 接收机钟差 | NaN | NaN | 整秒时 | 儒略日（整天） | 儒略日（小数天） 
    global GPST_UTC

    f1 = 1575.42 * 1E6;             % Hz
    f2 = 1227.60 * 1E6;             % Hz
    c = 2.99792458E+8;           % m/s
    lam1 = c / f1;
    lam2 = c / f2;
    phase1 = obsdata.L1 * lam1;
    phase2 = obsdata.L2 * lam2;
    phase = (f1^2 * phase1 - f2^2 * phase2) / (f1^2 - f2^2);
    p1 = obsdata.P1;
    p2 = obsdata.P2;
    p = obsdata.P;
    sp3_sat = sp3data.pos;
    clk_sat = clkdata.clk;
    number_of_sat = size(phase, 2);
    number_of_epo = size(phase, 1);
    
    phase(isnan(detect_pos)) = NaN;
    p1(isnan(detect_pos)) = NaN;
    p2(isnan(detect_pos)) = NaN;

    % 剔除卫星数少于4的历元
    for i = 1: number_of_epo
        num_NaN = numel(find(isnan(phase(i,:))));
        if number_of_sat - num_NaN < 4
            phase(i, ~isnan(phase(i,:))) = NaN;
        end
        if number_of_sat - numel(find(isnan(p1(i,:)))) < 4
            p1(i, ~isnan(p1(i,:))) = NaN;
        end
        if number_of_sat - numel(find(isnan(p2(i,:)))) < 4
            p2(i, ~isnan(p2(i,:))) = NaN;
        end
    end
    
    % 接收机钟差为NaN的置零
    spp_clk(isnan(spp_clk)) = 0;

    phasedata = NaN(number_of_sat*number_of_epo, 14);
    pseudodata = NaN(number_of_sat*number_of_epo, 14);
    t0 = 1;
    t1 = 1;
    for j = 1: number_of_sat
        for i = 1: number_of_epo
            if ~isnan(phase(i,j)) && ~isnan(detect_pos(i,j))
                JD1 = obsdata.epoch(i,1) + 2400000.5;
                JD2 = (obsdata.epoch(i,2) - GPST_UTC) / 86400;
                phasedata(t0, :) = [phase(i,j), obsdata.epochnum(i), detect_pos(i,j), j, sp3_sat(i,1:3,j),...
                    clk_sat(i,j), spp_clk(i), NaN, NaN, obsdata.epoch(i,3), JD1, JD2];
                t0 = t0 + 1;
                pseudodata(t1, :) = [p(i,j), obsdata.epochnum(i), NaN, j, sp3_sat(i,1:3,j),...
                    clk_sat(i,j), spp_clk(i), NaN, NaN, obsdata.epoch(i,3), JD1, JD2];
                t1 = t1 + 1;
            end
        end       
    end
    phasedata(isnan(phasedata(:,9)),:)=[];
    pseudodata(isnan(pseudodata(:,8)),:)=[];
end