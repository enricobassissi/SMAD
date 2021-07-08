%% --------------------------------------------------------------------- %% 
%% --------------------------- MARGIN ---------------------------------- %%
%% --------------------------------------------------------------------- %% 
% specific margins for LT
% MAR-DV-120 10% delta-v margin
% MAR-DV-130 Direct escape launch: An allocation of 45 m/s
% MAR-DV-150 interplanetary approach navigation manoeuvres 20 m/s
%% result after DT
clear; clc;
str_path=split(pwd, 'TrajOptimisation\direct transcript\follower');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
str_path_1=split(pwd, 'follower');
imp_path=string(str_path_1(1))+'functions';
addpath(genpath(imp_path));

%% 2RL final 15 mN
load('MA5_160kg_dry.mat')
% fuel_wet = 1-(exp(-sol.dV_tot*1e3/(g0*Isp))); 
g0 = 9.81;
Isp = 2600;

%% --- SC1
mass_depleted_Leg1 = MA.SC1.leg1.mass_depleted;
dV_associated_Leg1 = -g0*Isp*log(MA.SC1.leg1.mass_end/MA.SC1.leg1.mass_start); % -ve*ln(m_final/m_initial)
mass_depleted_Leg2 = MA.SC1.leg2.mass_depleted;
dV_associated_Leg2 = -g0*Isp*log(MA.SC1.leg2.mass_end/MA.SC1.leg2.mass_start); % -ve*ln(m_final/m_initial)
dV_not_margined_SC1 = [dV_associated_Leg1, dV_associated_Leg2];

% dV_margin_DV_010 =  dV_not_margined*0.05;
dV_margin_DV_010 =  [0,0];
% dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_080 = 0;
dV_margin_DV_120 = dV_not_margined_SC1*0.1;
dV_margin_DV_130 = 0.045; % launcher, direct injection
dV_margin_DV_150 = 0.020; % interplanetary approach
dV_margined_SC1 = [dV_not_margined_SC1(1)+dV_margin_DV_010(1)+dV_margin_DV_080+...
    dV_margin_DV_120(1)+dV_margin_DV_130,...
    dV_not_margined_SC1(2)+dV_margin_DV_010(2)+dV_margin_DV_120(2)+2*dV_margin_DV_150]
dV_margined_tot_SC1 = sum(dV_margined_SC1)
mass_frac_margined_SC1 = 1-(exp(-dV_margined_tot_SC1/(g0*Isp)))
mass_depleted_Leg1_margined = (1-(exp(-dV_margined_SC1(1)/(g0*Isp))))*MA.SC1.leg1.mass_start
mass_depleted_Leg2_margined = (1-(exp(-dV_margined_SC1(2)/(g0*Isp))))*MA.SC1.leg2.mass_start
mass_fuel_margined_SC1 = mass_depleted_Leg1_margined+mass_depleted_Leg2_margined

%% --- SC2
mass_depleted_Lega = MA.SC2.lega.mass_depleted;
dV_associated_Lega = -g0*Isp*log(MA.SC2.lega.mass_end/MA.SC2.lega.mass_start); % -ve*ln(m_final/m_initial)
mass_depleted_Legb = MA.SC2.legb.mass_depleted;
dV_associated_Legb = -g0*Isp*log(MA.SC2.legb.mass_end/MA.SC2.legb.mass_start); % -ve*ln(m_final/m_initial)
dV_not_margined_SC2 = [dV_associated_Lega, dV_associated_Legb];

% dV_margin_DV_010 =  dV_not_margined*0.05;
dV_margin_DV_010 =  [0,0];
% dV_margin_DV_080 = 0.030; % km/s % launcher
dV_margin_DV_080 = 0;
dV_margin_DV_120 = dV_not_margined_SC2*0.1;
dV_margin_DV_130 = 0.045; % launcher, direct injection
dV_margin_DV_150 = 0.020; % interplanetary approach
dV_margined_SC2 = [dV_not_margined_SC2(1)+dV_margin_DV_010(1)+dV_margin_DV_080+...
    dV_margin_DV_120(1)+dV_margin_DV_130,...
    dV_not_margined_SC2(2)+dV_margin_DV_010(2)+dV_margin_DV_120(2)+2*dV_margin_DV_150]
dV_margined_tot_SC2 = sum(dV_margined_SC2)
mass_frac_margined_SC2 = 1-(exp(-dV_margined_tot_SC2/(g0*Isp)))

mass_depleted_Lega_margined = (1-(exp(-dV_margined_SC2(1)/(g0*Isp))))*MA.SC2.lega.mass_start
mass_depleted_Legb_margined = (1-(exp(-dV_margined_SC2(2)/(g0*Isp))))*MA.SC2.legb.mass_start
mass_fuel_margined_SC2 = mass_depleted_Lega_margined+mass_depleted_Legb_margined
