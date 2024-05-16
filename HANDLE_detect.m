function [detect_pos, JUMP, pc_lc] = HANDLE_detect(obsdata)
% DESCRIPTION:     Detect.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    c = 2.99792458E+8;   % m/s
    f1 = 1575.42;   % MHz
    f2 = 1227.60;   % MHz
    
    lam1 = c / (f1 * 1E6);  % m
    lam2 = c / (f2 * 1E6);  % m
    
    L1 = obsdata.L1;
    L2 = obsdata.L2;
    P1 = obsdata.P1;
    P2 = obsdata.P2;
    P = obsdata.P;
    
    L1(L1==0) = NaN;
    L2(L2==0) = NaN;
    P1(P1==0) = NaN;
    P2(P2==0) = NaN;
    
    L = (f1^2 * L1 * lam1 - f2^2 * L2 * lam2) / (f1^2 - f2^2);
    pc_lc = P - L;
    
    % 宽巷组合
    MW = (L1 - L2) - (f1 - f2) / (f1 + f2) * (P1 / lam1 + P2 / lam2);
    
    phasedif = L - L2 * lam2;
    
    line = size(P1, 1);
    colu = size(P1, 2);
    detect_pos = NaN(line, colu);
    JUMP = zeros(line, colu);
    
    countpam = 0;
    aver_NWk_1 = 0;
    sigma2 = 0;
    aver_pclc = 0;
    for i = 1:colu
        k = 0;
        for j = 2: line
            if ~isnan(MW(j, i)) && ~isnan(MW(j-1, i))    
                
                k = k + 1;
                
                % 1.MW 组合
                aver_NWk = aver_NWk_1 + (MW(j, i) - aver_NWk_1) / k;
                sigma2_new = sigma2 + ((aver_NWk - aver_NWk_1)^2 - sigma2) / k;
                jump = MW(j, i) - aver_NWk_1;
                sigma4 = 4 * sqrt(sigma2_new);
                aver_NWk_1 = aver_NWk;
                sigma2 = sigma2_new;
                
                % 2.PC-LC
                aver_pclc = aver_pclc + (pc_lc(j,i) - aver_pclc) / k;
                
                % 认为发生周跳
                cs1 = j == 2;
                cs2 = abs(jump) > sigma4; 
                cs3 = abs(phasedif(j,i) - phasedif(j-1,i)) > 0.3; % ΔL1-L2
                cs4 = abs((pc_lc(j,i) - aver_pclc) / k) > 4;
                cs2 = 0;
                
                if cs1 || cs2 || cs3 || cs4
                    countpam = countpam + 1;
                    detect_pos(j, i) = countpam;
                    if cs1
                        detect_pos(j-1, i) = countpam;
                    end
                    
                    k = 1;
                    aver_NWk_1 = 0;
                    sigma2 = 0;
                    aver_pclc = 0;
                    if cs1
                        JUMP(j, i) = JUMP(j, i) + 1000;
                    end
                    if cs2
                        JUMP(j, i) = JUMP(j, i) + 200;
                    end
                    if cs3
                        JUMP(j, i) = JUMP(j, i) + 30;
                    end
                    if cs4
                        JUMP(j, i) = JUMP(j, i) + 4;
                    end
                else
                    detect_pos(j, i) = countpam;
                end
            elseif ~isnan(MW(j, i)) && isnan(MW(j-1, i))
                countpam = countpam + 1;
                detect_pos(j, i) = countpam;
                
                k = 1;
                aver_NWk_1 = MW(j,i);
                sigma2 = 0;
                aver_pclc = pc_lc(j,i);
            end
        end
    end
    
    % 周跳探测结果筛选：连续数据少于某阈值的段不处理
%     threshold = 10;
%     unique_nums = unique(detect_pos(~isnan(detect_pos)));
%     for i = 1:length(unique_nums)
%         target_num = unique_nums(i);
%         occurrences = sum(detect_pos(:) == target_num, 'omitnan');
%         if occurrences < threshold
%             detect_pos(detect_pos == target_num) = NaN;
%         end
%     end
    threshold = 10;
    [unique_nums, ~, id] = unique(detect_pos(~isnan(detect_pos)));
    occurrences = accumarray(id, 1);
    frequent = unique_nums(occurrences >= threshold);
    mask = ismember(detect_pos, frequent);
    detect_pos(~mask) = NaN;
end