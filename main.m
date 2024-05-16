%% 全局设置
HANDLE_gval;

%% 读取文件
global atxdata
% sp3data = HANDLE_readsp3('./input/20240103.SP3');
% clkdata = HANDLE_readclk("./input/20240103.CLK");
% obsdata = HANDLE_readobs("./input/0104-4006.rnx");
% atxdata = HANDLE_readatx("./input/igs14_2247.atx");
% dcbdata = HANDLE_readosb("./input/9.BIA");

sp3data = HANDLE_readsp3('./input/grace/igs22673.sp3');
clkdata = HANDLE_readclk("./input/grace/igs22673.clk_30s");
obsdata = HANDLE_readobs("./input/grace/GRCC1720.23O");
atxdata = HANDLE_readatx("./input/igs14_2247.atx");
dcbdata = HANDLE_readosb("./input/grace/COD0MGXFIN_20231720000_01D_01D_OSB.BIA");

%% 预处理
[sp3data, clkdata, obsdata] = HANDLE_interp(sp3data, clkdata, obsdata);
obsdata = HANDLE_combine(obsdata, dcbdata);
[detect_pos, JUMP, pc_lc] = HANDLE_detect(obsdata);

%% 单点定位
[orbit, obsdata] = spp(obsdata, clkdata, detect_pos);

%% 初轨确定
[init_orbit, JD2_init, orbits] = IOD(orbit);

%% 精密定轨
[phasedata, pseudodata] = preprocess_phase(obsdata, sp3data, clkdata, detect_pos, orbit(:,6), orbits(1,end));
[final_orbit, final_clk, final_amb, orbit_list, N_LOOP, X_correct, RMSE, inter_t, y_error] = phasepod(phasedata, pseudodata, init_orbit(1:6), JD2_init, orbits(1,end)-obsdata.epoch(1,3), obsdata.epoch(2,3)-obsdata.epoch(1,3));

%% 
orbits_gfz = readgfz('GNI1B_2023-06-21_C_04.txt');
yc_gfz = orbits_gfz(1:85011,2:4) - orbit_list(1:85011,1:3);
RTNerror = cmpRTN(yc_gfz, orbit_list(1:85011,1:6));
plot(RTNerror(:,1:end),'.')
legend('R','T','N')
mean(RTNerror)