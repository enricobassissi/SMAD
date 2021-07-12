%% --------------------------------------------------------------------- %%
%% ------------------ DT THE FLEXIBILITY ANALYSIS RESULT --------------- %%
%% --------------------------------------------------------------------- %%
%% Setup for default options
set(0, 'DefaultTextFontSize', 20)
set(0, 'DefaultAxesFontSize', 20)
set(0, 'DefaultLegendFontSize', 20)
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLineLineWidth', 1.8)
format short

%% Initializing the Environment
clear; close all; clc;

% Palette ESA
colors = [0    50   71;... % (1) DEEP SPACE
          207  29   57;... % (2) EXCITE RED +1
          0    142  122;... % (3) PURE TEAL +1
          251  171  24;... % (4) ENLIGHT YELLOW
          244  121  32;... % (5) ENLIGHT YELLOW +1
          150  1    54;... % (6) EXCITE RED +2
          167  85   52;... % (7) ENLIGHT YELLOW +2
          0    97   158;... % (8) TRUSTY AZURE +1
          30   51   120;... % (9) TRUSTY AZURE +2
          0    103  98;... % (10) PURE TEAL +2
          51   94   111;... % (11) DEEP SPACE -1
          0    174  157;... % (12) PURE TEAL
          0    0    0]./255; % (13) BLACK

%% add path of functions and python stuff
str_path=split(pwd, 'PostAnalysis\Flexibility');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
fl_path=split(pwd, 'Flexibility');
imp_path=string(fl_path(1));
addpath(genpath(imp_path));
dt_path=string(str_path(1))+'TrajOptimisation\direct transcript\follower';
addpath(genpath(dt_path));
dtf_path=string(str_path(1))+'TrajOptimisation\direct transcript\functions';
addpath(genpath(dtf_path));

%% Call to NASA JPL Horizons to get Asteroid's Ephemerides
% Import module of Python
try 
    module = py.importlib.import_module('neo_api_function');
catch
    copyfile(py_path+'\neo_api_function.py', pwd, 'f'); 
    module = py.importlib.import_module('neo_api_function');
end

%% simulation parameters
sim.mu_dim    = 132712440018              ; % actractor parameter [km^3 s^-2]
sim.DU        = 149597870.7               ; % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu        = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol     = 100;                    % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 

sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = -1;                     % direction of integration (1 FW, -1 BW), 
                                       % 1 is like imposing wet mass at beginning
sim.TOF_imposed_flag = 1;
sim.PS.Isp = 3200/sim.TU;  % non-dimensional specific impulse
% sim.PS.Isp = 4500/sim.TU;  % non-dimensional specific impulse % simone
sim.M1_end = 160; % SC wet mass [kg] %%
sim.M2_end = 160; % SC wet mass [kg] %%
sim.M_pods = 3.5; % mass of the payloads + landing stuff [kg] %%
sim.max_Available_Thrust = 0.02; % 5 [mN], BepiColombo is 250 mN but it's much bigger

%% %% LOAD INIT GUESS DATA AND FLEX DATA
load('160kg_dry_64mN.mat')
% load('ws_2RL_all_indietro_moo2.mat')
% the process there is at -1, from dry to wet... so we define the wet a posteriori
sol.M_start_SC1_leg1=output.m_SC1(1);
% sol.M_start_SC1_leg2=output.m_SC1(2*sim.n_sol);
sol.M_start_SC2_lega=output.m_SC2(1);
% sol.M_start_SC2_legb=output.m_SC2(2*sim.n_sol);
sol.Href = [output.Href.leg1,output.Href.leg2,output.Href.lega,output.Href.legb];
EphData = data;

clearvars -except colors sim sol r_encounter v_encounter EphData

load('LW_300_30.mat')
sim.direction = 1;

% clearvars -except colors sim sol r_encounter v_encounter EphData data 
%% ------------------------- Run Direct -------------------------------- %%
idx = 0;
for i = 1:length(lw_vector)
    idx = idx+1;
    fprintf('Cycle %d \n',idx);

    
    %% given that solution, build the NLI input variables
    MJD01 = flex{i}.MJD0*sim.TU/86400;
    % 1st spacecraft characteristic times
    MJDA1 = MJD01 + flex{i}.TOF1*sim.TU/86400; % mjd2000 passage of 1st sc on ast 1
    MJDD1 = MJDA1 + flex{i}.CT1*sim.TU/86400;
    MJDA2 = MJDD1 + flex{i}.TOF2*sim.TU/86400; % mjd2000 passage of 1st sc on ast 2

    % 2nd spacecraft characteristic times
    MJDAa = MJD01 + flex{i}.TOFa*sim.TU/86400; 
    MJDDa = MJDAa + flex{i}.CTa*sim.TU/86400;
    MJDAb = MJDDa + flex{i}.TOFb*sim.TU/86400; 

    % N REV1
    N_rev1 = data.Nrev1;
    % N REV2
    N_rev2 = data.Nrev2;
    % N REVa
    N_reva = data.Nreva;
    % N REVb
    N_revb = data.Nrevb;

    % C3 launcher
    v_inf_magn = flex{i}.v_inf_magn;
    az = flex{i}.az;
    elev = flex{i}.el;

    % asteroids
    asteroid_1 = "2020VV";
    asteroid_2 = "2009TD17";

    asteroid_a = "2011BP40";
    asteroid_b = "2021JE1";

    %% Computing position and velocity of the planets in that days
    % Departure from Earth
    [kep_EA,ksun] = uplanet(MJD01, 3);
    [r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
    r_EA = r_EA/sim.DU;
    v_EA = v_EA/sim.DU*sim.TU;

    % ARRIVAL at 1st ast
    [kep_ast_A1] = uNEO3(MJDA1,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
    [rA1, vA1] = sv_from_coe(kep_ast_A1,ksun); % km, km/s
    rA1 = rA1/sim.DU;
    vA1 = vA1/sim.DU*sim.TU;
    % DEPARTURE at 1st ast
    [kep_ast_D1] = uNEO3(MJDD1,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
    [rD1, vD1] = sv_from_coe(kep_ast_D1,ksun); % km, km/s
    rD1 = rD1/sim.DU;
    vD1 = vD1/sim.DU*sim.TU;

    % ARRIVAL at 2nd ast
    [kep_ast_A2] = uNEO3(MJDA2,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
    [rA2, vA2] = sv_from_coe(kep_ast_A2,ksun); % km, km/s
    rA2 = rA2/sim.DU;
    vA2 = vA2/sim.DU*sim.TU;

    % ARRIVAL at a_th ast
    [kep_ast_Aa] = uNEO3(MJDAa,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
    [rAa, vAa] = sv_from_coe(kep_ast_Aa,ksun); % km, km/s
    rAa = rAa/sim.DU;
    vAa = vAa/sim.DU*sim.TU;
    % DEPARTURE at a_th ast
    [kep_ast_Da] = uNEO3(MJDDa,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
    [rDa, vDa] = sv_from_coe(kep_ast_Da,ksun); % km, km/s
    rDa = rDa/sim.DU;
    vDa = vDa/sim.DU*sim.TU;

    % passage at b_th ast
    [kep_ast_Ab] = uNEO3(MJDAb,asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
    [rAb, vAb] = sv_from_coe(kep_ast_Ab,ksun); % km, km/s
    rAb = rAb/sim.DU;
    vAb = vAb/sim.DU*sim.TU;

    %% Launcher departure variable
    v_launcher = v_inf_magn*[cos(elev)*cos(az); cos(elev)*sin(az); sin(elev)];
    v_dep = v_EA + v_launcher/sim.DU*sim.TU;  %if parabolic escape (v_extra = 0)

    %% DT - SC 1 - LEG 1
    data.Tmax = 0.02; % max thrust available for EGI
    data.Pmax = 550; % max power related to 20 mN thrust on qitetiq t5
    data.Isp = 2900/sim.TU; % for 20 mN the Isp of Qinetiq T5 is 2900s
    data.ThrustModulationFlag = 1; % 0 -> nothing, only T_max
                                   % 1 -> only sun distance
                                   % 2 -> sun distance + time degradation
                                   % 3 -> sun distance + time degradation + SAA degradation
    data.n_int = 20; % discretisation selected
    data.angle_inplane_panels_max = 20; % degrees
    data.muS = sim.mu_dim;  
    data.MU = sol.mass_dry_and_pods_SC1+sol.tot_mass_depleted_SC1; % mass adimensionalisation on wet mass
    SC1.asteroid_1 = sol.asteroid_1;
    data.fmincon_iter = 100;

    SC1.leg1.v_in = v_dep;
    SC1.leg1.v_end = vA1;
    SC1.leg1.r_in = r_EA;
    SC1.leg1.r_end = rA1;
    SC1.leg1.N_rev = N_rev1;
    SC1.leg1.TOF = flex{i}.TOF1;

    DT_input = SC1.leg1;
    [SC1.leg1.HS, SC1.leg1.RES, SC1.leg1.timead, SC1.leg1.Xad, SC1.leg1.Xpropad, ...
        SC1.leg1.Xpropad_DT] = DT_executable_flexibility(DT_input, sim, data);

    % -- new mass after DT calculations
    SC1.leg1.Href = SC1.leg1.RES.Href;
    SC1.leg1.mass_start = SC1.leg1.HS.X(1,7);
    SC1.leg1.mass_end = SC1.leg1.HS.X(end,7);
    SC1.leg1.mass_depleted = SC1.leg1.mass_start - SC1.leg1.mass_end;
    SC1.leg1.mass_fraction = SC1.leg1.mass_depleted/SC1.leg1.mass_start;

    SC1.leg1.PROPULSION = nikita(SC1.leg1, data, sim, colors);
    SC1.leg1.EPS = marco_elia(SC1.leg1,r_encounter.EA,data,sim,SC1.leg1.Href,colors);

    flex_dt{i}.SC1.leg1 = SC1.leg1;

    %% DT - SC 1 - LEG 2
    data.Tmax = 0.015; % max thrust available for EGI
    data.Pmax = 450; % max power related to 15 mN thrust on qitetiq t5
    data.Isp = 2750/sim.TU; % for 15 mN the Isp of Qinetiq T5 is 2750s
    data.ThrustModulationFlag = 2;
    data.angle_inplane_panels_max = 10; % degrees
    data.n_int = 20;
    % data.MU = sol.M_start_SC1_leg2; % mass adimensionalisation on wet mass
    data.MU = flex_dt{i}.SC1.leg1.mass_end - sim.M_pods; % mass adimensionalisation on wet mass
    SC1.asteroid_2 = sol.asteroid_2;

    SC1.leg2.v_in = vD1;
    SC1.leg2.v_end = vA2;
    SC1.leg2.r_in = rD1;
    SC1.leg2.r_end = rA2;
    SC1.leg2.N_rev = N_rev2;
    SC1.leg2.TOF = flex{i}.TOF2;

    DT_input = SC1.leg2;
    [SC1.leg2.HS, SC1.leg2.RES, SC1.leg2.timead, SC1.leg2.Xad, SC1.leg2.Xpropad, ...
        SC1.leg2.Xpropad_DT] = DT_executable_flexibility(DT_input, sim, data);

    % -- new mass after DT calculations
    SC1.leg2.Href = SC1.leg2.RES.Href;
    SC1.leg2.mass_start = SC1.leg2.HS.X(1,7);
    SC1.leg2.mass_end = SC1.leg2.HS.X(end,7);
    SC1.leg2.mass_depleted = SC1.leg2.mass_start - SC1.leg2.mass_end;
    SC1.leg2.mass_fraction = SC1.leg2.mass_depleted/SC1.leg2.mass_start;

    SC1.leg2.PROPULSION = nikita(SC1.leg2, data, sim, colors);
    SC1.leg2.EPS = marco_elia(SC1.leg2,r_encounter.astD1,data,sim,SC1.leg2.Href,colors);

    flex_dt{i}.SC1.leg2 = SC1.leg2;
    
    %% DT - SC 2 - LEG a
    data.Tmax = 0.015; % max thrust available for EGI
    data.Pmax = 450; % max power related to 15 mN thrust on qitetiq t5
    data.Isp = 2750/sim.TU; % for 15 mN the Isp of Qinetiq T5 is 2750s
    data.ThrustModulationFlag = 2;
    data.angle_inplane_panels_max = 10; % degrees
    data.n_int = 20;
    data.MU = sol.mass_dry_and_pods_SC2+sol.tot_mass_depleted_SC2; % mass adimensionalisation on wet mass
    SC2.asteroid_a = sol.asteroid_a;

    SC2.lega.v_in = v_dep;
    SC2.lega.v_end = vAa;
    SC2.lega.r_in = r_EA;
    SC2.lega.r_end = rAa;
    SC2.lega.N_rev = N_reva;
    SC2.lega.TOF = flex{i}.TOFa;

    DT_input = SC2.lega;
    [SC2.lega.HS, SC2.lega.RES, SC2.lega.timead, SC2.lega.Xad, SC2.lega.Xpropad, ...
        SC2.lega.Xpropad_DT] = DT_executable_flexibility(DT_input, sim, data);

    % -- new mass after DT calculations
    SC2.lega.Href = SC2.lega.RES.Href;
    SC2.lega.mass_start = SC2.lega.HS.X(1,7);
    SC2.lega.mass_end = SC2.lega.HS.X(end,7);
    SC2.lega.mass_depleted = SC2.lega.mass_start - SC2.lega.mass_end;
    SC2.lega.mass_fraction = SC2.lega.mass_depleted/SC2.lega.mass_start;

    SC2.lega.PROPULSION = nikita(SC2.lega, data, sim, colors);
    SC2.lega.EPS = marco_elia(SC2.lega,r_encounter.EA,data,sim,SC2.lega.Href,colors);

    flex_dt{i}.SC2.lega = SC2.lega;
    
    %% DT - SC 2 - LEG b
    data.Tmax = 0.015; % max thrust available for EGI
    data.Pmax = 450; % max power related to 15 mN thrust on qitetiq t5
    data.Isp = 2750/sim.TU; % for 15 mN the Isp of Qinetiq T5 is 2750s
    data.ThrustModulationFlag = 1;
    data.angle_inplane_panels_max = 20; % degrees
    data.n_int = 20;
    % data.MU = sol.M_start_SC2_legb; % mass adimensionalisation on wet mass
    data.MU = flex_dt{i}.SC2.lega.mass_end - sim.M_pods; % mass adimensionalisation on wet mass
    SC2.asteroid_b = sol.asteroid_b;

    SC2.legb.v_in = vDa;
    SC2.legb.v_end = vAb;
    SC2.legb.r_in = rDa;
    SC2.legb.r_end = rAb;
    SC2.legb.N_rev = N_revb;
    SC2.legb.TOF = flex{i}.TOFb;

    DT_input = SC2.legb;
    [SC2.legb.HS, SC2.legb.RES, SC2.legb.timead, SC2.legb.Xad, SC2.legb.Xpropad, ...
        SC2.legb.Xpropad_DT] = DT_executable_flexibility(DT_input, sim, data);

    % -- new mass after DT calculations
    SC2.legb.Href = SC2.legb.RES.Href;
    SC2.legb.mass_start = SC2.legb.HS.X(1,7);
    SC2.legb.mass_end = SC2.legb.HS.X(end,7);
    SC2.legb.mass_depleted = SC2.legb.mass_start - SC2.legb.mass_end;
    SC2.legb.mass_fraction = SC2.legb.mass_depleted/SC2.legb.mass_start;

    SC2.legb.PROPULSION = nikita(SC2.legb, data, sim, colors);
    SC2.legb.EPS = marco_elia(SC2.legb,r_encounter.astDa,data,sim,SC2.legb.Href,colors);
    
    flex_dt{i}.SC2.legb = SC2.legb;
    
    close all
    
end
