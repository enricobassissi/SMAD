%% --------------------------------------------------------------------- %%
%% ----- GIVEN THE SOLUTION OF A WS OBTAINED BY GENETIC ALGORITHM ------ %%
%% ---- THE RESULT IS PASSED TO THE DIRECT TRASCRIPTION TO OPTIMISE ---- %%
%% ----------- THE THRUST DISTRIBUTION AND THE MASS CONSUMPTION -------- %%
%% --------------------------------------------------------------------- %%
%% ------- ALL THE LEGS ARE SINGOLARLY PASSED THROUGH THE FUNCTION ----- %%
%% --------------------------------------------------------------------- %%
%% Default options
clear; close all; clc;

set(0, 'DefaultTextFontSize', 20) % modify it if too small
set(0, 'DefaultAxesFontSize', 20) % modify it if too small
set(0, 'DefaultLegendFontSize', 20) % modify it if too small
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLineLineWidth', 1.8)

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
          0    0    0]./255; % (12) BLACK

%% work environment setup
str_path=split(pwd, 'TrajOptimisation\direct transcript\follower');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
str_path_1=split(pwd, 'follower');
imp_path=string(str_path_1(1))+'functions';
addpath(genpath(imp_path));

%% LOAD INIT GUESS DATA
load('160kg_dry_64mN.mat')
% load('ws_2RL_all_indietro_moo2.mat')
% the process there is at -1, from dry to wet... so we define the wet a posteriori
sol.M_start_SC1_leg1=output.m_SC1(1);
% sol.M_start_SC1_leg2=output.m_SC1(2*sim.n_sol);
sol.M_start_SC2_lega=output.m_SC2(1);
% sol.M_start_SC2_legb=output.m_SC2(2*sim.n_sol);
sim.direction = 1;
sol.Href = [output.Href.leg1,output.Href.leg2,output.Href.lega,output.Href.legb];
EphData = data;

clearvars -except colors sim sol r_encounter v_encounter EphData

%% --------------------- DIRECT TRANSCRIPTION -------------------------- %%
%% Thruster
% --- 13 mN max
% data.Isp = 2600/sim.TU; % for 13 mN the Isp of nikita thruster is 2600s
% --- 15 mN max
% data.Tmax = 0.015; % max thrust available for EGI
% data.Pmax = 450; % max power related to 15 mN thrust on qitetiq t5
% data.Isp = 2750/sim.TU; % for 15 mN the Isp of Qinetiq T5 is 2750s
% --- 20 mN max
data.Tmax = 0.02; % max thrust available for EGI
data.Pmax = 550; % max power related to 20 mN thrust on qitetiq t5
data.Isp = 2900/sim.TU; % for 20 mN the Isp of Qinetiq T5 is 2900s
% --- 25 mN max
% data.Isp = 3200/sim.TU; % for 25 mN the Isp of RitEvo T5 is 3200s
data.LifeTime = 20*1e3; % h, for Qinetiq T5
data.TotImpulse_max = 3.5*1e6; % [Ns], Qinetiq T5
% --- Hall Thrusers
% data.Tmax = 0.025; % max thrust available for Hall thruster
% data.Isp = 1600/sim.TU; % for 25 mN the Isp of SETS ST40 Hall Thruster is 1600s
% data.LifeTime = 5000; % h, for hall thruster
% data.TotImpulse_max = data.Tmax*data.LifeTime*3600; % [Ns], Hall thruster
data.muS = sim.mu_dim;   
data.ThrustModulationFlag = 1; % 0 -> nothing, only T_max
                               % 1 -> only sun distance
                               % 2 -> sun distance + time degradation
                               % 3 -> sun distance + time degradation + SAA degradation

%% DT - SC 1 - LEG 1
data.Tmax = 0.02; % max thrust available for EGI
data.Pmax = 550; % max power related to 20 mN thrust on qitetiq t5
data.Isp = 2900/sim.TU; % for 20 mN the Isp of Qinetiq T5 is 2900s
data.n_int = 20; % discretisation selected
data.angle_inplane_panels_max = 20; % degrees
data.MU = sol.M_start_SC1_leg1; % mass adimensionalisation on wet mass
SC1.asteroid_1 = sol.asteroid_1;

v_launcher = sol.v_inf_magn/sim.DU*sim.TU*[cos(sol.el)*cos(sol.az); cos(sol.el)*sin(sol.az); sin(sol.el)];
v_dep = v_encounter.EA + v_launcher;  %if parabolic escape (v_extra = 0)
SC1.leg1.v_in = v_dep;
SC1.leg1.v_end = v_encounter.astA1;
SC1.leg1.r_in = r_encounter.EA;
SC1.leg1.r_end = r_encounter.astA1;
SC1.leg1.N_rev = sol.Nrev(1);
SC1.leg1.TOF = sol.TOF1_ADIM;
SC1.leg1.Href = sol.Href(:,1);

DT_input = SC1.leg1;
tic
[SC1.leg1.HS, SC1.leg1.RES, SC1.leg1.timead, SC1.leg1.Xad, SC1.leg1.Xpropad, ...
    SC1.leg1.Xpropad_DT] = DT_executable(DT_input, sim, data);
el_time.DT_SC1_leg1 = toc/60; % minutes

SC1.leg1.PROPULSION = nikita(SC1.leg1, data, sim, colors);
SC1.leg1.EPS = marco_elia(SC1.leg1,r_encounter.EA,data,sim,SC1.leg1.Href,colors);

%% -- if run with less than 100 points
data.n_int = 100;
[SC1.leg1] = interpoliamo_porc(SC1.leg1,data.n_int);
SC1.leg1.v_in = v_dep;
SC1.leg1.v_end = v_encounter.astA1;
SC1.leg1.r_in = r_encounter.EA;
SC1.leg1.r_end = r_encounter.astA1;
SC1.leg1.N_rev = sol.Nrev(1);
SC1.leg1.TOF = sol.TOF1_ADIM;
SC1.leg1.Href = sol.Href(:,1);
SC1.mjd2000_dep = sol.departure_mjd2000;
SC1.TOF1 = sol.TOF1; 
SC1.TOF2 = sol.TOF2;
SC1.CT1 = sol.CT1;

%% analysis on the results
IDcase = 'SC1 - LEG1';
DT_plots(SC1.leg1.HS, SC1.leg1.RES, SC1.leg1.timead, SC1.leg1.Xad, SC1.leg1.Xpropad, ...
    SC1.leg1.Xpropad_DT, IDcase, colors);

% -- new mass after DT calculations
SC1.leg1.mass_start = SC1.leg1.HS.X(1,7);
SC1.leg1.mass_end = SC1.leg1.HS.X(end,7);
SC1.leg1.mass_depleted = SC1.leg1.mass_start - SC1.leg1.mass_end;
SC1.leg1.mass_fraction = SC1.leg1.mass_depleted/SC1.leg1.mass_start;

SC1.leg1.PROPULSION = nikita(SC1.leg1, data, sim, colors);
SC1.leg1.EPS = marco_elia(SC1.leg1,r_encounter.EA,data,sim,SC1.leg1.Href,colors);

% % -- coasting between the two legs
[SC1.coasting.leg1] = relative_position_coasting_stuff(SC1.asteroid_1,SC1.mjd2000_dep,...
    SC1.TOF1,SC1.CT1,data.n_int,EphData,sim,colors);%leg1.timead(end)*sim.TU/86400

%% DT - SC 1 - LEG 2
data.Tmax = 0.015; % max thrust available for EGI
data.Pmax = 450; % max power related to 15 mN thrust on qitetiq t5
data.Isp = 2750/sim.TU; % for 15 mN the Isp of Qinetiq T5 is 2750s
data.ThrustModulationFlag = 2;
data.angle_inplane_panels_max = 10; % degrees
data.n_int = 70;
% data.MU = sol.M_start_SC1_leg2; % mass adimensionalisation on wet mass
data.MU = SC1.leg1.mass_end - sim.M_pods; % mass adimensionalisation on wet mass
SC1.asteroid_2 = sol.asteroid_2;

SC1.leg2.v_in = v_encounter.astD1;
SC1.leg2.v_end = v_encounter.astA2;
SC1.leg2.r_in = r_encounter.astD1;
SC1.leg2.r_end = r_encounter.astA2;
SC1.leg2.N_rev = sol.Nrev(2);
SC1.leg2.TOF = sol.TOF2_ADIM;
SC1.leg2.Href = sol.Href(:,2);

DT_input = SC1.leg2;
tic
[SC1.leg2.HS, SC1.leg2.RES, SC1.leg2.timead, SC1.leg2.Xad, SC1.leg2.Xpropad, ...
    SC1.leg2.Xpropad_DT] = DT_executable(DT_input, sim, data);
el_time.DT_SC1_leg2 = toc/60; % minutes

SC1.leg2.PROPULSION = nikita(SC1.leg2, data, sim, colors);
SC1.leg2.EPS = marco_elia(SC1.leg2,r_encounter.astD1,data,sim,SC1.leg2.Href,colors);

%% -- if run with less than 100 points
data.n_int = 100;
[SC1.leg2] = interpoliamo_porc(SC1.leg2,data.n_int);
SC1.leg2.v_in = v_encounter.astD1;
SC1.leg2.v_end = v_encounter.astA2;
SC1.leg2.r_in = r_encounter.astD1;
SC1.leg2.r_end = r_encounter.astA2;
SC1.leg2.N_rev = sol.Nrev(2);
SC1.leg2.TOF = sol.TOF2_ADIM;
SC1.leg2.Href = sol.Href(:,2);
SC1.mjd2000_dep = sol.departure_mjd2000;
SC1.TOF1 = sol.TOF1; 
SC1.TOF2 = sol.TOF2;
SC1.CT1 = sol.CT1;

%% analysis on the results
IDcase = 'SC1 - LEG2';
DT_plots(SC1.leg2.HS, SC1.leg2.RES, SC1.leg2.timead, SC1.leg2.Xad, SC1.leg2.Xpropad, ...
    SC1.leg2.Xpropad_DT, IDcase, colors);

SC1.leg2.PROPULSION = nikita(SC1.leg2, data, sim, colors);
SC1.leg2.EPS = marco_elia(SC1.leg2,r_encounter.astD1,data,sim,SC1.leg2.Href,colors);

% -- new mass after DT calculations
SC1.leg2.mass_start = SC1.leg2.HS.X(1,7);
SC1.leg2.mass_end = SC1.leg2.HS.X(end,7);
SC1.leg2.mass_depleted = SC1.leg2.mass_start - SC1.leg2.mass_end;
SC1.leg2.mass_fraction = SC1.leg2.mass_depleted/SC1.leg2.mass_start;

% -- coasting on the last asteroid, same coasting time
SC1.mjd2000_dep_ast1 = SC1.mjd2000_dep+SC1.TOF1+SC1.CT1;
[SC1.coasting.leg2] = relative_position_coasting_stuff(SC1.asteroid_2,SC1.mjd2000_dep_ast1,...
    SC1.TOF2,SC1.CT1,data.n_int,EphData,sim,colors); %SC1.leg2.timead(end)*sim.TU/86400

%% -- building of a uniform position, time, thrust vector for SC1
SC1.uniform.R = [SC1.leg1.EPS.R_cartesian; SC1.coasting.leg1.r_ast; 
    SC1.leg2.EPS.R_cartesian; SC1.coasting.leg2.r_ast;]./sim.DU; % AU
figure()
% plot3(SC1.uniform.R(:,1),SC1.uniform.R(:,2),SC1.uniform.R(:,3));
plot3(SC1.uniform.R(1:data.n_int,1),SC1.uniform.R(1:data.n_int,2),SC1.uniform.R(1:data.n_int,3),...
    'Color',colors(1,:),'DisplayName','Leg 1');
hold on
plot3(SC1.uniform.R(data.n_int+1:2*data.n_int,1),SC1.uniform.R(data.n_int+1:2*data.n_int,2),...
    SC1.uniform.R(data.n_int+1:2*data.n_int,3),...
    'Color',colors(2,:),'DisplayName','Rendezvous 1');
plot3(SC1.uniform.R(2*data.n_int+1:3*data.n_int,1),SC1.uniform.R(2*data.n_int+1:3*data.n_int,2),...
    SC1.uniform.R(2*data.n_int+1:3*data.n_int,3),...
    'Color',colors(3,:),'DisplayName','Leg 2');
plot3(SC1.uniform.R(3*data.n_int+1:4*data.n_int,1),SC1.uniform.R(3*data.n_int+1:4*data.n_int,2),...
    SC1.uniform.R(3*data.n_int+1:4*data.n_int,3),...
    'Color',colors(4,:),'DisplayName','Rendezvous 2');
plot3(r_encounter.EA(1), r_encounter.EA(2), r_encounter.EA(3),'*',...
    'Color',colors(8,:),'DisplayName','EA Dep');
plot3(r_encounter.astA1(1), r_encounter.astA1(2), r_encounter.astA1(3),'^',...
    'DisplayName','Ast 1 Arr');
plot3(r_encounter.astD1(1), r_encounter.astD1(2), r_encounter.astD1(3),'^',...
    'DisplayName','Ast 1 Dep');
plot3(r_encounter.astA2(1), r_encounter.astA2(2), r_encounter.astA2(3),'^',...
    'DisplayName','Ast 2 Arr');
plot3(SC1.coasting.leg2.r_ast(end,1)./sim.DU, SC1.coasting.leg2.r_ast(end,2)./sim.DU,...
    SC1.coasting.leg2.r_ast(end,3)./sim.DU,'^', 'DisplayName','End Mission');
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); axis equal, grid on;
legend('show');

SC1.uniform.time = [SC1.mjd2000_dep+SC1.leg1.timead.*sim.TU/86400; SC1.coasting.leg1.time;
    SC1.coasting.leg1.time(end)+SC1.leg2.timead.*sim.TU/86400; SC1.coasting.leg2.time;]; % mjd 2000

SC1.uniform.Thrust = [SC1.leg1.EPS.T_transf_orbit; zeros(data.n_int,3);
    SC1.leg2.EPS.T_transf_orbit; zeros(data.n_int,3)]; % N
figure()
plot(SC1.uniform.time, SC1.uniform.Thrust)
xlabel('t [mjd2000]'); ylabel('T [N]'); legend('Tx','Ty','Tz');
figure()
plot(SC1.uniform.time, vecnorm(SC1.uniform.Thrust,2,2))
xlabel('t [mjd2000]'); ylabel('T magnitude [N]');

%% DT - SC 2 - LEG a
data.Tmax = 0.015; % max thrust available for EGI
data.Pmax = 450; % max power related to 15 mN thrust on qitetiq t5
data.Isp = 2750/sim.TU; % for 15 mN the Isp of Qinetiq T5 is 2750s
data.ThrustModulationFlag = 2;
data.angle_inplane_panels_max = 10; % degrees
data.n_int = 70;
data.MU = sol.M_start_SC2_lega; % mass adimensionalisation on wet mass
SC2.asteroid_a = sol.asteroid_a;

SC2.lega.v_in = v_dep;
SC2.lega.v_end = v_encounter.astAa;
SC2.lega.r_in = r_encounter.EA;
SC2.lega.r_end = r_encounter.astAa;
SC2.lega.N_rev = sol.Nrev(3);
SC2.lega.TOF = sol.TOFa_ADIM;
SC2.lega.Href = sol.Href(:,3);

DT_input = SC2.lega;
tic
[SC2.lega.HS, SC2.lega.RES, SC2.lega.timead, SC2.lega.Xad, SC2.lega.Xpropad, ...
    SC2.lega.Xpropad_DT] = DT_executable(DT_input, sim, data);
el_time.DT_SC2_lega = toc/60; % minutes

SC2.lega.PROPULSION = nikita(SC2.lega, data, sim, colors);
SC2.lega.EPS = marco_elia(SC2.lega,r_encounter.EA,data,sim,SC2.lega.Href,colors);

%% -- if run with less than 100 points
data.n_int = 100;
[SC2.lega] = interpoliamo_porc(SC2.lega,data.n_int);
SC2.lega.v_in = v_dep;
SC2.lega.v_end = v_encounter.astAa;
SC2.lega.r_in = r_encounter.EA;
SC2.lega.r_end = r_encounter.astAa;
SC2.lega.N_rev = sol.Nrev(3);
SC2.lega.TOF = sol.TOFa_ADIM;
SC2.lega.Href = sol.Href(:,3);
SC2.mjd2000_dep = sol.departure_mjd2000;
SC2.TOFa = sol.TOFa; 
SC2.TOFb = sol.TOFb;
SC2.CTa = sol.CTa;

%% analysis on the results
IDcase = 'SC2 - LEGa';
DT_plots(SC2.lega.HS, SC2.lega.RES, SC2.lega.timead, SC2.lega.Xad, SC2.lega.Xpropad, ...
    SC2.lega.Xpropad_DT, IDcase, colors);

SC2.lega.PROPULSION = nikita(SC2.lega, data, sim, colors);
SC2.lega.EPS = marco_elia(SC2.lega,r_encounter.EA,data,sim,SC2.lega.Href,colors);

% -- new mass after DT calculations
SC2.lega.mass_start = SC2.lega.HS.X(1,7);
SC2.lega.mass_end = SC2.lega.HS.X(end,7);
SC2.lega.mass_depleted = SC2.lega.mass_start - SC2.lega.mass_end;
SC2.lega.mass_fraction = SC2.lega.mass_depleted/SC2.lega.mass_start;

% -- coasting between the two legs
[SC2.coasting.lega] = relative_position_coasting_stuff(SC2.asteroid_a,SC2.mjd2000_dep,...
    SC2.TOFa,SC2.CTa,data.n_int,EphData,sim,colors);

%% DT - SC 2 - LEG b
data.Tmax = 0.015; % max thrust available for EGI
data.Pmax = 450; % max power related to 15 mN thrust on qitetiq t5
data.Isp = 2750/sim.TU; % for 15 mN the Isp of Qinetiq T5 is 2750s
data.ThrustModulationFlag = 1;
data.angle_inplane_panels_max = 20; % degrees
data.n_int = 20;
% data.MU = sol.M_start_SC2_legb; % mass adimensionalisation on wet mass
data.MU = SC2.lega.mass_end - sim.M_pods; % mass adimensionalisation on wet mass
SC2.asteroid_b = sol.asteroid_b;

SC2.legb.v_in = v_encounter.astDa;
SC2.legb.v_end = v_encounter.astAb;
SC2.legb.r_in = r_encounter.astDa;
SC2.legb.r_end = r_encounter.astAb;
SC2.legb.N_rev = sol.Nrev(4);
SC2.legb.TOF = sol.TOFb_ADIM;
SC2.legb.Href = sol.Href(:,4);

DT_input = SC2.legb;
tic
[SC2.legb.HS, SC2.legb.RES, SC2.legb.timead, SC2.legb.Xad, SC2.legb.Xpropad, ...
    SC2.legb.Xpropad_DT] = DT_executable(DT_input, sim, data);
el_time.DT_SC2_legb = toc/60; % minutes
el_time.DT_tot = el_time.DT_SC1_leg1+el_time.DT_SC1_leg2+el_time.DT_SC2_lega+el_time.DT_SC2_legb; % min

SC2.legb.PROPULSION = nikita(SC2.legb, data, sim, colors);
SC2.legb.EPS = marco_elia(SC2.legb,r_encounter.astDa,data,sim,SC2.legb.Href,colors);

%% -- if run with less than 100 points
data.n_int = 100;
[SC2.legb] = interpoliamo_porc(SC2.legb,data.n_int);
SC2.legb.v_in = v_encounter.astDa;
SC2.legb.v_end = v_encounter.astAb;
SC2.legb.r_in = r_encounter.astDa;
SC2.legb.r_end = r_encounter.astAb;
SC2.legb.N_rev = sol.Nrev(4);
SC2.legb.TOF = sol.TOFb_ADIM;
SC2.legb.Href = sol.Href(:,4);
SC2.mjd2000_dep = sol.departure_mjd2000;
SC2.TOFa = sol.TOFa; 
SC2.TOFb = sol.TOFb;
SC2.CTa = sol.CTa;

%% post analysis on the results
IDcase = 'SC2 - LEGb';
DT_plots(SC2.legb.HS, SC2.legb.RES, SC2.legb.timead, SC2.legb.Xad, SC2.legb.Xpropad, ...
    SC2.legb.Xpropad_DT, IDcase, colors);

SC2.legb.PROPULSION = nikita(SC2.legb, data, sim, colors);
SC2.legb.EPS = marco_elia(SC2.legb,r_encounter.astDa,data,sim,SC2.legb.Href,colors);

% -- new mass after DT calculations
SC2.legb.mass_start = SC2.legb.HS.X(1,7);
SC2.legb.mass_end = SC2.legb.HS.X(end,7);
SC2.legb.mass_depleted = SC2.legb.mass_start - SC2.legb.mass_end;
SC2.legb.mass_fraction = SC2.legb.mass_depleted/SC2.legb.mass_start;

% -- coasting on the last asteroid, same coasting time
SC2.mjd2000_dep_asta = SC2.mjd2000_dep+SC2.TOFa+SC2.CTa;
[SC2.coasting.legb] = relative_position_coasting_stuff(SC2.asteroid_b,SC2.mjd2000_dep_asta,...
    SC2.TOFb,SC2.CTa,data.n_int,EphData,sim,colors);

%% -- building of a uniform position, time, thrust vector for SC2
SC2.uniform.R = [SC2.lega.EPS.R_cartesian; SC2.coasting.lega.r_ast;
    SC2.legb.EPS.R_cartesian; SC2.coasting.legb.r_ast]./sim.DU; % AU
figure()
% plot3(SC2.uniform.R(:,1),SC2.uniform.R(:,2),SC2.uniform.R(:,3));
plot3(SC2.uniform.R(1:data.n_int,1),SC2.uniform.R(1:data.n_int,2),SC2.uniform.R(1:data.n_int,3),...
    'Color',colors(1,:),'DisplayName','Leg 1');
hold on
plot3(SC2.uniform.R(data.n_int+1:2*data.n_int,1),SC2.uniform.R(data.n_int+1:2*data.n_int,2),...
    SC2.uniform.R(data.n_int+1:2*data.n_int,3),...
    'Color',colors(2,:),'DisplayName','Rendezvous 1');
plot3(SC2.uniform.R(2*data.n_int+1:3*data.n_int,1),SC2.uniform.R(2*data.n_int+1:3*data.n_int,2),...
    SC2.uniform.R(2*data.n_int+1:3*data.n_int,3),...
    'Color',colors(3,:),'DisplayName','Leg 2');
plot3(SC2.uniform.R(3*data.n_int+1:4*data.n_int,1),SC2.uniform.R(3*data.n_int+1:4*data.n_int,2),...
    SC2.uniform.R(3*data.n_int+1:4*data.n_int,3),...
    'Color',colors(4,:),'DisplayName','Rendezvous 2');
plot3(r_encounter.EA(1), r_encounter.EA(2), r_encounter.EA(3),'*',...
    'Color',colors(8,:),'DisplayName','EA Dep');
plot3(r_encounter.astAa(1), r_encounter.astAa(2), r_encounter.astAa(3),'^',...
    'DisplayName','Ast a Arr');
plot3(r_encounter.astDa(1), r_encounter.astDa(2), r_encounter.astDa(3),'^',...
    'DisplayName','Ast a Dep');
plot3(r_encounter.astAb(1), r_encounter.astAb(2), r_encounter.astAb(3),'^',...
    'DisplayName','Ast 2 Arr');
plot3(SC2.coasting.legb.r_ast(end,1)./sim.DU, SC2.coasting.legb.r_ast(end,2)./sim.DU,...
    SC2.coasting.legb.r_ast(end,3)./sim.DU,'^', 'DisplayName','End Mission');
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); axis equal, grid on;
legend('show');

SC2.uniform.time = [SC2.mjd2000_dep+SC2.lega.timead.*sim.TU/86400; SC2.coasting.lega.time;
    SC2.coasting.lega.time(end)+SC2.legb.timead.*sim.TU/86400; SC2.coasting.legb.time]; % [mjd 2000]

SC2.uniform.Thrust = [SC2.lega.EPS.T_transf_orbit; zeros(data.n_int,3);
    SC2.legb.EPS.T_transf_orbit; zeros(data.n_int,3)]; % [N]
figure()
plot(SC2.uniform.time, SC2.uniform.Thrust)
xlabel('t [mjd2000]'); ylabel('T [N]'); legend('Tx','Ty','Tz');
figure()
plot(SC2.uniform.time, vecnorm(SC2.uniform.Thrust,2,2))
xlabel('t [mjd2000]'); ylabel('T magnitude [N]');

%% --------------------------------------------------------------------- %%
%% ----------------- MISSION ANALYSIS OUTPUT STRUCT -------------------- %%
%% --------------------------------------------------------------------- %%
SC1_temp = SC1;
SC2_temp = SC2;
SC1_temp.leg1 = rmfield(SC1_temp.leg1, 'RES');
SC1_temp.leg1 = rmfield(SC1_temp.leg1, 'Xad');
SC1_temp.leg1 = rmfield(SC1_temp.leg1, 'Xpropad');
SC1_temp.leg1 = rmfield(SC1_temp.leg1, 'Xpropad_DT');
SC1_temp.leg2 = rmfield(SC1_temp.leg2, 'RES');
SC1_temp.leg2 = rmfield(SC1_temp.leg2, 'Xad');
SC1_temp.leg2 = rmfield(SC1_temp.leg2, 'Xpropad');
SC1_temp.leg2 = rmfield(SC1_temp.leg2, 'Xpropad_DT');

SC2_temp.lega = rmfield(SC2_temp.lega, 'RES');
SC2_temp.lega = rmfield(SC2_temp.lega, 'Xad');
SC2_temp.lega = rmfield(SC2_temp.lega, 'Xpropad');
SC2_temp.lega = rmfield(SC2_temp.lega, 'Xpropad_DT');
SC2_temp.legb = rmfield(SC2_temp.legb, 'RES');
SC2_temp.legb = rmfield(SC2_temp.legb, 'Xad');
SC2_temp.legb = rmfield(SC2_temp.legb, 'Xpropad');
SC2_temp.legb = rmfield(SC2_temp.legb, 'Xpropad_DT');

MA.SC1 = SC1_temp;
MA.SC2 = SC2_temp;

%% dimensionalisation
% -- sc1
MA.SC1.leg1.timead = MA.SC1.leg1.timead*sim.TU;
MA.SC1.leg1.EPS.time = MA.SC1.leg1.EPS.time*sim.TU;
MA.SC1.leg1.TOF = MA.SC1.leg1.TOF*sim.TU/86400;
MA.SC1.coasting.leg1.time = (MA.SC1.coasting.leg1.time - (MA.SC1.mjd2000_dep+MA.SC1.TOF1))*86400;

MA.SC1.leg2.timead = MA.SC1.leg2.timead*sim.TU;
MA.SC1.leg2.EPS.time = MA.SC1.leg2.EPS.time*sim.TU;
MA.SC1.leg2.TOF = MA.SC1.leg2.TOF*sim.TU/86400;
MA.SC1.coasting.leg2.time = (MA.SC1.coasting.leg2.time - (MA.SC1.mjd2000_dep+MA.SC1.TOF1+MA.SC1.CT1+MA.SC1.TOF2))*86400;

MA.SC1.uniform.R = MA.SC1.uniform.R.*sim.DU;
MA.SC1.uniform.time = (MA.SC1.uniform.time - MA.SC1.mjd2000_dep)*86400;

% -- sc2
MA.SC2.lega.timead = MA.SC2.lega.timead*sim.TU;
MA.SC2.lega.EPS.time = MA.SC2.lega.EPS.time*sim.TU;
MA.SC2.lega.TOF = MA.SC2.lega.TOF*sim.TU/86400;
MA.SC2.coasting.lega.time = (MA.SC2.coasting.lega.time - (MA.SC2.mjd2000_dep+MA.SC2.TOFa))*86400;

MA.SC2.legb.timead = MA.SC2.legb.timead*sim.TU;
MA.SC2.legb.EPS.time = MA.SC2.legb.EPS.time*sim.TU;
MA.SC2.legb.TOF = MA.SC2.legb.TOF*sim.TU/86400;
MA.SC2.coasting.legb.time = (MA.SC2.coasting.legb.time - (MA.SC2.mjd2000_dep+MA.SC2.TOFa+MA.SC2.CTa+MA.SC2.TOFb))*86400;

MA.SC2.uniform.R = MA.SC2.uniform.R.*sim.DU;
MA.SC2.uniform.time = (MA.SC2.uniform.time - MA.SC2.mjd2000_dep)*86400;

%% obj encounter
r_encounter_dim.EA = r_encounter.EA*sim.DU;
r_encounter_dim.astA1 = r_encounter.astA1*sim.DU;
r_encounter_dim.astD1 = r_encounter.astD1*sim.DU;
r_encounter_dim.astA2 = r_encounter.astA2*sim.DU;
r_encounter_dim.astAa = r_encounter.astAa*sim.DU;
r_encounter_dim.astDa = r_encounter.astDa*sim.DU;
r_encounter_dim.astAb = r_encounter.astAb*sim.DU;
MA.r_encounter = r_encounter_dim;

%% plot with MA struct
MA.data.n_int = data.n_int;
figure()
% plot3(MA.SC2.uniform.R(:,1),MA.SC2.uniform.R(:,2),MA.SC2.uniform.R(:,3));
plot3(MA.SC2.uniform.R(1:data.n_int,1),MA.SC2.uniform.R(1:data.n_int,2),MA.SC2.uniform.R(1:data.n_int,3),...
    'Color',colors(1,:),'DisplayName','Leg 1');
hold on
plot3(MA.SC2.uniform.R(data.n_int+1:2*data.n_int,1),MA.SC2.uniform.R(data.n_int+1:2*data.n_int,2),...
    MA.SC2.uniform.R(data.n_int+1:2*data.n_int,3),...
    'Color',colors(2,:),'DisplayName','Rendezvous 1');
plot3(MA.SC2.uniform.R(2*data.n_int+1:3*data.n_int,1),MA.SC2.uniform.R(2*data.n_int+1:3*data.n_int,2),...
    MA.SC2.uniform.R(2*data.n_int+1:3*data.n_int,3),...
    'Color',colors(3,:),'DisplayName','Leg 2');
plot3(MA.SC2.uniform.R(3*data.n_int+1:4*data.n_int,1),MA.SC2.uniform.R(3*data.n_int+1:4*data.n_int,2),...
    MA.SC2.uniform.R(3*data.n_int+1:4*data.n_int,3),...
    'Color',colors(4,:),'DisplayName','Rendezvous 2');
plot3(MA.r_encounter.EA(1), MA.r_encounter.EA(2), MA.r_encounter.EA(3),'*',...
    'Color',colors(8,:),'DisplayName','EA Dep');
plot3(MA.r_encounter.astAa(1), MA.r_encounter.astAa(2), MA.r_encounter.astAa(3),'^',...
    'DisplayName','Ast a Arr');
plot3(MA.r_encounter.astDa(1), MA.r_encounter.astDa(2), MA.r_encounter.astDa(3),'^',...
    'DisplayName','Ast a Dep');
plot3(MA.r_encounter.astAb(1), MA.r_encounter.astAb(2), MA.r_encounter.astAb(3),'^',...
    'DisplayName','Ast 2 Arr');
plot3(MA.SC2.coasting.legb.r_ast(end,1), MA.SC2.coasting.legb.r_ast(end,2),...
    MA.SC2.coasting.legb.r_ast(end,3),'^', 'DisplayName','End Mission');
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); axis equal, grid on;
legend('show');

figure()
plot(MA.SC2.uniform.time, MA.SC2.uniform.Thrust)
xlabel('t [s]'); ylabel('T [N]'); legend('Tx','Ty','Tz');
figure()
plot(MA.SC2.uniform.time, vecnorm(MA.SC2.uniform.Thrust,2,2))
xlabel('t [s]'); ylabel('T magnitude [N]');

