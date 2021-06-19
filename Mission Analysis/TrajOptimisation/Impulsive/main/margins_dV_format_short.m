%% --------------------------------------------------------------------- %% 
%% --------------------------- MARGIN ---------------------------------- %%
%% --------------------------------------------------------------------- %% 
%% add path of functions and python stuff
str_path=split(pwd, 'TrajOptimisation\Impulsive\main');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'main');
imp_path=string(str_path(1));
addpath(genpath(imp_path));

%% 1+4 FLYBY
clear; clc;
% load('soo_ps_flyby_3-33_1dayeachpointTraj.mat');
load('ws_1FI_4.3_vrel3.mat');
dV_not_margined = [sol.dV_extra_launch,sol.dVast1,sol.dVast2,sol.dVast3];
check_if_bigger_than_10ms = dV_not_margined*0.05; % any case
dV_margin_DV_010 = sum(check_if_bigger_than_10ms);
dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_090 = 0.015; % km/s % GA
dV_margined = sum(dV_not_margined)+dV_margin_DV_010+dV_margin_DV_080


%% 1 SC
clear; clc;
load('soo_ps_rendezvous_18-51.mat');
dV_not_margined = [sol.dV_extra_launch,sol.dV2,sol.dV3,sol.dV4,sol.dV5,sol.dV6,sol.dV7,sol.dV8];
check_if_bigger_than_10ms = dV_not_margined*0.05; % any case
dV_margin_DV_010 = sum(check_if_bigger_than_10ms);
dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_090 = 0.015; % km/s % GA
dV_margined = sum(dV_not_margined)+dV_margin_DV_010+dV_margin_DV_080

%% 2SC - FLYBY
clear; clc;
load('soo_ps_flyby_2sc_5-14.mat');
dV_not_margined = [sol.dV_extra_launch_sc1,sol.dV_extra_launch_sc2,sol.dVast1,sol.dVasta];
check_if_bigger_than_10ms = dV_not_margined*0.05; % any case
dV_margin_DV_010 = sum(check_if_bigger_than_10ms);
dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_090 = 0.015; % km/s % GA
dV_margined = sum(dV_not_margined)+dV_margin_DV_010+dV_margin_DV_080

%% 2SC - RV
clear; clc;
% load('soo_ps_rendezvous_2sc_23-97.mat');
load('ws_2RL_10.mat');
dV_not_margined = [sol.dV1_ast12, sol.dV2_ast12, sol.dV1_astab, sol.dV2_astab, sol.dV2_EAasta_sc2, sol.dV2_EAast1_sc1, sol.dV_extra_launch_sc2, sol.dV_extra_launch_sc1];
check_if_bigger_than_10ms = dV_not_margined*0.05; % any case
dV_margin_DV_010 = sum(check_if_bigger_than_10ms);
dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_090 = 0.015; % km/s % GA
dV_margined = sum(dV_not_margined)+dV_margin_DV_010+dV_margin_DV_080

%% 1+4 - GA EA
clear; clc;
load('soo_ps_flyby_GA_4-74.mat');
dV_not_margined = [sol.delta_V_p, sol.dVast1, sol.dVast2, sol.dVast3];
check_if_bigger_than_10ms = dV_not_margined*0.05; % any case
dV_margin_DV_010 = sum(check_if_bigger_than_10ms);
dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_090 = 0.015; % km/s % GA
dV_margined = sum(dV_not_margined)+dV_margin_DV_010+dV_margin_DV_080+dV_margin_DV_090

%% 1+4 - GA MA
clear; clc;
load('soo_ps_flyby_MARS_GA_11-85.mat');
dV_not_margined = [sol.dVast3, sol.dVast2, sol.dVast1, sol.delta_V_p, sol.dV_extra_launch];
check_if_bigger_than_10ms = dV_not_margined*0.05; % any case
dV_margin_DV_010 = sum(check_if_bigger_than_10ms);
dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_090 = 0.015; % km/s % GA
dV_margined = sum(dV_not_margined)+dV_margin_DV_010+dV_margin_DV_080+dV_margin_DV_090

%% 2SC - FLYBY - 0.24
clear; clc;
load('soo_ps_flyby_2sc_0-24.mat');
dV_not_margined = [sol.dVasta, sol.dVast1];
check_if_bigger_than_10ms = dV_not_margined*0.05; % any case
dV_margin_DV_010 = sum(check_if_bigger_than_10ms);
dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_090 = 0.015; % km/s % GA
dV_margined = sum(dV_not_margined)+dV_margin_DV_010+dV_margin_DV_080
