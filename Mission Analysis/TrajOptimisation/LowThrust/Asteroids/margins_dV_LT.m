%% --------------------------------------------------------------------- %% 
%% --------------------------- MARGIN ---------------------------------- %%
%% --------------------------------------------------------------------- %% 
% specific margins for LT
% MAR-DV-120 10% delta-v margin
% MAR-DV-130 Direct escape launch: An allocation of 45 m/s
% MAR-DV-150 interplanetary approach navigation manoeuvres 20 m/s

%% 1FL
clear; clc;
str_path=split(pwd, 'TrajOptimisation\LowThrust\Asteroids');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'Asteroids');
imp_path=string(str_path(1));
addpath(genpath(imp_path));

load(imp_path+'Asteroids\Workspaces\ws_1FL_0.6_pods_45mN.mat')
% fuel_wet = 1-(exp(-sol.dV_tot*1e3/(g0*Isp))); 
g0 = 9.81;
Isp = 3200;
dV_not_margined = [sol.dV_associated_Leg1, sol.dV_associated_Leg2, sol.dV_associated_Leg3, sol.dV_associated_Leg4];
% dV_margin_DV_010 =  dV_not_margined*0.05;
dV_margin_DV_010 =  [0,0,0,0];
% dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_080 = 0;
dV_margin_DV_120 = dV_not_margined*0.1;
dV_margin_DV_130 = 0.045; % launcher, direct injection
dV_margin_DV_150 = 0.020; % interplanetary approach
dV_margined = [dV_not_margined(1)+dV_margin_DV_010(1)+dV_margin_DV_080+...
    dV_margin_DV_120(1)+dV_margin_DV_130,...
    dV_not_margined(2)+dV_margin_DV_010(2)+dV_margin_DV_120(2)+dV_margin_DV_150,...
    dV_not_margined(3)+dV_margin_DV_010(3)+dV_margin_DV_120(3)+dV_margin_DV_150,...
    dV_not_margined(4)+dV_margin_DV_010(4)+dV_margin_DV_120(4)+dV_margin_DV_150]

mass_frac_margined = 1-(exp(-sum(dV_margined)/(g0*Isp)))

%% 1RL
clear; clc;
str_path=split(pwd, 'TrajOptimisation\LowThrust\Asteroids');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'Asteroids');
imp_path=string(str_path(1));
addpath(genpath(imp_path));

load(imp_path+'Asteroids\Workspaces\ws_1RL_0.7_27mN.mat')
% fuel_wet = 1-(exp(-sol.dV_tot*1e3/(g0*Isp))); 
g0 = 9.81;
Isp = 3200;
dV_not_margined = [sol.dV_associated_Leg1, sol.dV_associated_Leg2, sol.dV_associated_Leg3, sol.dV_associated_Leg4];
% dV_margin_DV_010 =  dV_not_margined*0.05;
dV_margin_DV_010 =  [0,0,0,0];
% dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_080 = 0;
dV_margin_DV_120 = dV_not_margined*0.1;
dV_margin_DV_130 = 0.045; % launcher, direct injection
dV_margin_DV_150 = 0.020; % interplanetary approach
dV_margined = [dV_not_margined(1)+dV_margin_DV_010(1)+dV_margin_DV_080+...
    dV_margin_DV_120(1)+dV_margin_DV_130,...
    dV_not_margined(2)+dV_margin_DV_010(2)+dV_margin_DV_120(2)+dV_margin_DV_150,...
    dV_not_margined(3)+dV_margin_DV_010(3)+dV_margin_DV_120(3)+dV_margin_DV_150,...
    dV_not_margined(4)+dV_margin_DV_010(4)+dV_margin_DV_120(4)+dV_margin_DV_150]
sum(dV_margined)
mass_frac_margined = 1-(exp(-sum(dV_margined)/(g0*Isp)))

%% 2FL
clear; clc;
str_path=split(pwd, 'TrajOptimisation\LowThrust\Asteroids');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'Asteroids');
imp_path=string(str_path(1));
addpath(genpath(imp_path));

load(imp_path+'Asteroids\Workspaces\ws_2FL_0.3_pods_16mN.mat')
% fuel_wet = 1-(exp(-sol.dV_tot*1e3/(g0*Isp))); 
g0 = 9.81;
Isp = 3200;
dV_not_margined_SC1 = [sol.dV_associated_Leg1, sol.dV_associated_Leg2]
% dV_margin_DV_010 =  dV_not_margined*0.05;
dV_margin_DV_010 =  [0,0];
% dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_080 = 0;
dV_margin_DV_120 = dV_not_margined_SC1*0.1;
dV_margin_DV_130 = 0.045; % launcher, direct injection
dV_margin_DV_150 = 0.020; % interplanetary approach
dV_margined_SC1 = [dV_not_margined_SC1(1)+dV_margin_DV_010(1)+dV_margin_DV_080+...
    dV_margin_DV_120(1)+dV_margin_DV_130,...
    dV_not_margined_SC1(2)+dV_margin_DV_010(2)+dV_margin_DV_120(2)+dV_margin_DV_150]
sum(dV_margined_SC1)
mass_frac_margined_SC1 = 1-(exp(-sum(dV_margined_SC1)/(g0*Isp)))

dV_not_margined_SC2 = [sol.dV_associated_Lega, sol.dV_associated_Legb]
% dV_margin_DV_010 =  dV_not_margined*0.05;
dV_margin_DV_010 =  [0,0];
% dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_080 = 0;
dV_margin_DV_120 = dV_not_margined_SC2*0.1;
dV_margin_DV_130 = 0.045; % launcher, direct injection
dV_margin_DV_150 = 0.020; % interplanetary approach
dV_margined_SC2 = [dV_not_margined_SC2(1)+dV_margin_DV_010(1)+dV_margin_DV_080+...
    dV_margin_DV_120(1)+dV_margin_DV_130,...
    dV_not_margined_SC2(2)+dV_margin_DV_010(2)+dV_margin_DV_120(2)+dV_margin_DV_150]
sum(dV_margined_SC2)
mass_frac_margined_SC2 = 1-(exp(-sum(dV_margined_SC2)/(g0*Isp)))

%% 2RL
clear; clc;
str_path=split(pwd, 'TrajOptimisation\LowThrust\Asteroids');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'Asteroids');
imp_path=string(str_path(1));
addpath(genpath(imp_path));

load(imp_path+'Asteroids\Ale Workspaces\case22_for_margins.mat')
% fuel_wet = 1-(exp(-sol.dV_tot*1e3/(g0*Isp))); 
g0 = 9.81;
Isp = 3200;
dV_not_margined_SC1 = [sol.dV_associated_Leg1, sol.dV_associated_Leg2+1000]
% dV_margin_DV_010 =  dV_not_margined*0.05;
dV_margin_DV_010 =  [0,0];
% dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_080 = 0;
dV_margin_DV_120 = dV_not_margined_SC1*0.1;
dV_margin_DV_130 = 0.045; % launcher, direct injection
dV_margin_DV_150 = 0.020; % interplanetary approach
dV_margined_SC1 = [dV_not_margined_SC1(1)+dV_margin_DV_010(1)+dV_margin_DV_080+...
    dV_margin_DV_120(1)+dV_margin_DV_130,...
    dV_not_margined_SC1(2)+dV_margin_DV_010(2)+dV_margin_DV_120(2)+dV_margin_DV_150]
sum(dV_margined_SC1)
mass_frac_margined_SC1 = 1-(exp(-sum(dV_margined_SC1)/(g0*Isp)))

dV_not_margined_SC2 = [sol.dV_associated_Lega+1000, sol.dV_associated_Legb]
% dV_margin_DV_010 =  dV_not_margined*0.05;
dV_margin_DV_010 =  [0,0];
% dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_080 = 0;
dV_margin_DV_120 = dV_not_margined_SC2*0.1;
dV_margin_DV_130 = 0.045; % launcher, direct injection
dV_margin_DV_150 = 0.020; % interplanetary approach
dV_margined_SC2 = [dV_not_margined_SC2(1)+dV_margin_DV_010(1)+dV_margin_DV_080+...
    dV_margin_DV_120(1)+dV_margin_DV_130,...
    dV_not_margined_SC2(2)+dV_margin_DV_010(2)+dV_margin_DV_120(2)+dV_margin_DV_150]
sum(dV_margined_SC2)
mass_frac_margined_SC2 = 1-(exp(-sum(dV_margined_SC2)/(g0*Isp)))