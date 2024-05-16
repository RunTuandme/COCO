function [rd, r_gnss_ITRS_row, r_rho] = TOOL_lighttime(t0, JD1, JD2, activeS, r_leopre_GCRS, v_leopre_GCRS, clk_state, HG)
% DESCRIPTION:     Calculate light time.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           t0, second.
%                  JD1, the large part of the Julian Date.
%                  JD2, the small part of the Julian Date.
%                  activeS, Active GNSS satellites in search.
%                  r_leopre_GCRS, satellite position in GCRF.
%                  v_leopre_GCRS, satellite velocity in GCRF.
%                  clk_state, receiver clock offset.
%                  HG, rotation matrix from GCRF to ITRF.
    global GM
    rd = 0;
    c = 2.99792458E+8;
    leo_clk = clk_state / c;
    t0_gnss = t0 - leo_clk;
    while true      
        [r_gnss_ITRS_row, v_gnss_ITRS_row] = TOOL_getSP3posvel(rd, t0_gnss, JD1, JD2-rd/c/86400, activeS, HG);
%         distance = norm(r_gnss_ITRS_row - r_leopre_ITRS');
        eop = TOOL_geteop(JD1-2400000.5+JD2-rd/c/86400);
        [r_gnss_GCRS, v_gnss_GCRS, HGnew] = itrs2gcrs(r_gnss_ITRS_row', v_gnss_ITRS_row', JD1, JD2-rd/c/86400, eop.xp, eop.yp, eop.LOD, eop.UT1_UTC, eop.dX, eop.dY);
        HG = HGnew;
%         VG = [-r_leopre_ITRS(2); r_leopre_ITRS(1); 0] * omega_earth; % 地球自转修正
%         RG = r_gnss_ITRS_row' - r_leopre_ITRS;
%         r_leopre_GCRS = HG' * r_leopre_ITRS;
        distance = norm(r_gnss_GCRS - r_leopre_GCRS);
        d_rela = 2 * GM / c^3 * log((norm(r_leopre_GCRS) + norm(r_gnss_GCRS) + distance) / (norm(r_leopre_GCRS) + norm(r_gnss_GCRS) - distance)) + 2 * dot(r_gnss_GCRS, v_gnss_GCRS) / c^2 ...
            - 2 * dot(r_leopre_GCRS, v_leopre_GCRS) / c^2;
%         rdnew = distance - c * TOOL_getCLK(t0, activeS) + clk_state + d_rela * c - RG' * VG / c;
        rdnew = distance + d_rela * c;
        if abs(rd - rdnew) < 0.1
            break
        end
        rd = rdnew;
    end
    r_rho = r_leopre_GCRS - r_gnss_GCRS;
end

