function [orbit, obsdata] = spp(obsdata, clkdata, detect_pos)
% DESCRIPTION:     Single Point Positioning.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% OUTPUT:          orbit(format), JD1 | JD2 | X | Y | Z | clk | RMS | t0
    global GPST_UTC
    p = obsdata.P;
    clk = clkdata.clk;
    NE = size(p, 1);
    orbit = NaN(NE, 8);
%     POS_satfixed = NaN(NE, 3);
  
    for i = 1: NE
        if mod(i,1000)==0
            disp(['SPP解算至第 ' num2str(i) ' 个历元'])
        end
        JD1 = obsdata.epoch(i,1) + 2400000.5;
        JD2 = (obsdata.epoch(i,2)-GPST_UTC) / 86400;
        t0 = obsdata.epoch(i,3);
        line = p(i, :);
        validObservations = ~isnan(clk(i,:)) & ~isnan(line) & ~isnan(detect_pos(i,:));
        activeS = find(validObservations == 1);   
        P = line(validObservations);

        % Check for sufficient observations
        if size(P,2) < 4
            continue
        end
        C = clk(i, validObservations);
        
        try
        [pos, cut_activeS] = LST(P', C', t0, JD1, JD2, activeS); 
        catch
            a=0;
        end
        if isnan(pos(1))
            obsdata.P1(i,:) = NaN;
            obsdata.P2(i,:) = NaN;
            obsdata.P(i,:) = NaN;
            obsdata.L1(i,:) = NaN;
            obsdata.L2(i,:) = NaN;
            obsdata.L(i,:) = NaN;
            continue
        end
        if ~isempty(cut_activeS)
            obsdata.P1(i,cut_activeS) = NaN;
            obsdata.P2(i,cut_activeS) = NaN;
            obsdata.P(i,cut_activeS) = NaN;
            obsdata.L1(i,cut_activeS) = NaN;
            obsdata.L2(i,cut_activeS) = NaN;
            obsdata.L(i,cut_activeS) = NaN;
        end
        orbit(i, :) = [JD1 JD2 pos t0];
%         POS_satfixed(i, :) = [pps(1) pps(2) pps(3)];
    end
    dlmwrite('./output/spp.dat', orbit);
end

function [pos, cut_activeS] = LST(rho, clk, t0, JD1, JD2, activeS)
    % activeS： prn
    global GM
    omega_earth = 7.292115146706979e-5; % rad/s
    c = 2.99792458E+8;
    N = size(rho,1);
    S = ones(4, 1);
    HG = sofa_C2T(JD1, JD2);
    [POS, VEL] = TOOL_getSP3posvel(zeros(N,1), t0, JD1, JD2, activeS, HG);
    
    rms = Inf;
    rms_new = 1E6;
    counti = 0;
    cut_activeS = [];
    
    while (abs(rms_new - rms) > 0.01) || counti < 10      
        if counti >= 5
            [POS, VEL] = TOOL_getSP3posvel(r, t0, JD1, JD2, activeS, HG);
        end

        % GNSS PV
        X = POS(:,1);
        Y = POS(:,2);
        Z = POS(:,3);
        VX= VEL(:,1);
        VY= VEL(:,2);
        VZ= VEL(:,3);
        
        % Check elevation angle
        if counti >= 5
            [f3d, satfixed2earthfixed] = TOOL_getapcoffset(JD1, JD2, S(1:3), HG);
            S(1:3) = S(1:3) + f3d;
            cutgnss = [];
            for i = 1:length(activeS)
                POS_satfixed = satfixed2earthfixed' * [X(i)-S(1); Y(i)-S(2); Z(i)-S(3)];
                dxy = sqrt((POS_satfixed(2))^2 + (POS_satfixed(1))^2);
                if -POS_satfixed(3) / dxy < tan(5 * pi / 180)
                    cutgnss = [cutgnss; i];
                end
            end
%             if ~isempty(cutgnss)
%                 cut_activeS = [cut_activeS activeS(cutgnss)];
%                 rho(cutgnss) = [];
%                 activeS(cutgnss) = [];
%                 clk(cutgnss) = [];
%                 X(cutgnss) = [];
%                 Y(cutgnss) = [];
%                 Z(cutgnss) = [];
%                 VX(cutgnss) = [];
%                 VY(cutgnss) = [];
%                 VZ(cutgnss) = [];
%                 N = size(rho,1);
%             end
        end
        
        if isempty(rho) || size(rho, 1) < 4
            pos = NaN(1,5);
            return
        end
        
        VG = [-S(2); S(1); 0] * omega_earth; % Earth rotation correction
        RG = [X-S(1) Y-S(2) Z-S(3)];
        
        if counti < 5
            r = sqrt((X - S(1)).^2 + (Y - S(2)).^2 + (Z - S(3)).^2) - c * clk - RG * VG / c + S(4);
        else
            r_leo = HG' * [S(1); S(2); S(3)];
            d_rela = zeros(length(activeS), 1);
            for i = 1:length(activeS)
                eop = TOOL_geteop(JD1-2400000.5+JD2-r(i)/c/86400);
                r_ITRF = [X(i); Y(i); Z(i)];
                v_ITRF = [VX(i); VY(i); VZ(i)];
                [r_GCRS, v_GCRS] = itrs2gcrs(r_ITRF, v_ITRF, JD1, JD2-r(i)/c/86400, eop.xp, eop.yp, eop.LOD, eop.UT1_UTC, eop.dX, eop.dY);
                d_rela(i) = 2 * GM / c^3 * log((norm(r_leo + r_GCRS) + rho(i)) / (norm(r_leo + r_GCRS) - rho(i))) + ...
            2 * dot(r_GCRS, v_GCRS) / c^2;
            end
            clk = TOOL_getCLK(t0 - r / c, activeS);
            r = sqrt((X - S(1)).^2 + (Y - S(2)).^2 + (Z - S(3)).^2) - c * clk + d_rela * c - RG * VG / c + S(4);
        end
       
        lx = (X - S(1)) ./ r;
        ly = (Y - S(2)) ./ r;
        lz = (Z - S(3)) ./ r; 
        G = [-lx -ly -lz ones(N,1)];   
        b = rho - r;
        DX = (G' * G) \ (G' * b);
        rms = rms_new;
        rms_new = sqrt(sum(b.^2) / N);
        S = S + DX;
        counti = counti + 1;
        
        if counti > 9
            % Outlier removal
            sigma3 = 3 * sqrt(b'*b / length(b));
%             cutgnss = [];
%             for i = 1:length(b)
%                 if abs(b(i)) > sigma3 || abs(b(i)) > 10
%                     cutgnss = [cutgnss; i];
%                 end
%             end
            cutgnss = abs(b) > sigma3;
            
            if ~isempty(cutgnss) && any(cutgnss)
                rho(cutgnss) = [];
                activeS(cutgnss) = [];
                clk(cutgnss) = [];
                X(cutgnss) = [];
                Y(cutgnss) = [];
                Z(cutgnss) = [];
                VX(cutgnss) = [];
                VY(cutgnss) = [];
                VZ(cutgnss) = [];
                N = size(rho,1);
                rms_new = Inf;
                counti = 5;
            end
            if size(rho,1) < 4
                S(1:4) = NaN(1,4);
                rms_new = NaN;
                break
            end
        end
        
        if counti > 100
            break
        end
    end
    pos = [S(1) S(2) S(3) S(4) rms_new];
end