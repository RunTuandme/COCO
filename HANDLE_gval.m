function HANDLE_gval
% DESCRIPTION:     Define the golbal variables.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    % global JD1_pod       % Defined in function : IOD
    global obs_stepwithfile 
    obs_stepwithfile = 1;  % obs in use / obs from file

    global obscodes
    % 'G': GPS, 'C':BDS, 'E':GAL, 'GC':GPS&BDS
    obscodes.sys = 'G';
    obscodes.OG = {'P1', 'L1', 'P2', 'L2'};
%     obscodes.OG = {'C1C', 'L1C', 'C2W', 'L2P'};
    obscodes.OC = {'C2I', 'L2I', 'C6I', 'L6I'};
    obscodes.OE = {'C1X', 'L1P', 'C5I', 'L5P'};
    obscodes.G1 = 'C1W';
    obscodes.G2 = 'C2W';
    obscodes.Gnum = 32;
    obscodes.Cnum = 46;
    obscodes.Enum = 36;
    obscodes.num = obscodes.Gnum + obscodes.Cnum + obscodes.Enum;
    
    global exp_switch
    exp_switch = false;

    global xls
    global xpl
    global TERM
    cd ./SOFA
    load('xls.mat','-mat','xls');
    load('xpl.mat','-mat','xpl');
    load('TERM.mat','-mat','TERM');
    cd ..

    % gain from IERS, or set zero.
    % see: https://www.iers.org/IERS/EN/DataProducts/tools/eop_of_today/eop_of_today_tool.html
    global EOPmat
    global interxp
    global interyp
    global interUT1_TUC
    global interLOD
    global interdX
    global interdY
    EOPmat = HANDLE_readEOP('./EOP/finals.all.iau2000.txt');
    xprow = ~isnan(EOPmat(:, 2));
    yprow = ~isnan(EOPmat(:, 3));
    UT1_TUCrow = ~isnan(EOPmat(:, 4));
    LODrow = ~isnan(EOPmat(:, 5));
    dXrow = ~isnan(EOPmat(:, 6));
    dYrow = ~isnan(EOPmat(:, 7));
    interxp = griddedInterpolant(EOPmat(xprow,1), EOPmat(xprow, 2), 'spline');
    interyp = griddedInterpolant(EOPmat(yprow,1), EOPmat(yprow, 3), 'spline');
    interUT1_TUC = griddedInterpolant(EOPmat(UT1_TUCrow,1), EOPmat(UT1_TUCrow, 4), 'spline');
    interLOD = griddedInterpolant(EOPmat(LODrow,1), EOPmat(LODrow, 5), 'spline');
    interdX = griddedInterpolant(EOPmat(dXrow,1), EOPmat(dXrow, 6), 'spline');
    interdY = griddedInterpolant(EOPmat(dYrow,1), EOPmat(dYrow, 7), 'spline');
    
    global PC
    cd ./DE
    load DE440Coeff.mat DE440Coeff
    PC = DE440Coeff;
    cd ..
    
    global gravity_type      % 1: static EGM08 2: static EIGEN GL04C 3: GOCO06s
    gravity_type = 3;
    
    global clm0
    global clm
    global slm
    global fes2004
    global GM
    global ae
    if gravity_type == 1
        cd ./EGM2008
        load('EGM08.mat','-mat','GM','ae','clm0','clm','slm');
        cd ..
    elseif gravity_type == 2
        cd ./EIGENGL04C
        load('EIGENGL04.mat','-mat','GM','ae','clm0','clm','slm');
        cd ..
    elseif gravity_type == 3
        cd ./GOCO06s
        global tc0 tc ts acosc0 acosc acoss asinc0 asinc asins
        load('GOCO06s.mat','-mat','GM','ae','clm0','clm','slm','tc0','tc','ts','acosc0','acosc','acoss','asinc0','asinc','asins');
        cd ..
    end
        
    
    cd ./FES2004
    load('matlabfes2004.mat','-mat','fes2004');
    cd ..
    
    global GPST_UTC
    global TAI_UTC
    global TT_TAI
    
    GPST_UTC = 18;
    TAI_UTC = 37;
    TT_TAI = 32.184;
    
    global apcoffset 
    apcoffset = [0.2602;  -0.0013;  -0.4770];
end