%% --------------------------------------------------------------------- %%
%% ---- APPLY TO THE SOLUTION OF NOMINAL DIRECT AN ERROR OF ATTITUDE --- %%
%% ----- OF +/-1° ON THE IN PLANE AND OUT OF PLANE THRUST ANGLES ------- %%
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
str_path=split(pwd, 'TrajOptimisation\direct transcript\functions');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
str_path_1=split(pwd, 'functions');
imp_path=string(str_path_1(1))+'follower';
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
                               
%% --------------------- Alessio stuff -------------------- %%
%% DT - SC 1 - LEG 1
data.Tmax = 0.02; % max thrust available for EGI
data.Pmax = 550; % max power related to 20 mN thrust on qitetiq t5
data.Isp = 2900/sim.TU; % for 20 mN the Isp of Qinetiq T5 is 2900s
data.n_int = 20; % discretisation selected
data.angle_inplane_panels_max = 20; % degrees
data.MU = sol.M_start_SC1_leg1; % mass adimensionalisation on wet mass
SC1.asteroid_1 = sol.asteroid_1;
data.fmincon_iter = 500;

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
[SC1.leg1.HS, SC1.leg1.ANGLES_AS, SC1.leg1.timead, SC1.leg1.Xad, ...
    SC1.leg1.Xpropad_DT, SC1.leg1.Xpropad_AS] = DT_executable_attitude_sensitivity(DT_input, ...
    sim, data);

%% -- if run with less than 100 points
n = 100;
method = 'linear';
t_new = linspace(SC1.leg1.timead(1),SC1.leg1.timead(end),n)';
SC1.temp.timead       = t_new;
SC1.temp.HS.X = interp1(SC1.leg1.timead,SC1.leg1.HS.X,t_new,method);
SC1.temp.Xpropad_AS = interp1(SC1.leg1.timead,SC1.leg1.Xpropad_AS,t_new,method);
SC1.temp.Xpropad_DT = interp1(SC1.leg1.timead,SC1.leg1.Xpropad_DT,t_new,method);
SC1.temp.HS.beta = interp1(SC1.leg1.timead,SC1.leg1.HS.beta,t_new,method);
SC1.temp.HS.alpha = interp1(SC1.leg1.timead,SC1.leg1.HS.alpha,t_new,method);
SC1.temp.ANGLES_AS = interp1(SC1.leg1.timead,SC1.leg1.ANGLES_AS,t_new,method);
SC1.temp.mass_end = SC1.leg1.HS.X(end,7);

SC1.leg1 = SC1.temp;
SC1 = rmfield(SC1, 'temp');

% ---- cartesian orbit
% ---- nominal propagation
r_transf_orbit_nom  = [SC1.leg1.HS.X(:,1).*cos(SC1.leg1.HS.X(:,2)), ...
    SC1.leg1.HS.X(:,1).*sin(SC1.leg1.HS.X(:,2)), SC1.leg1.HS.X(:,3)];
R_cartesian_nom = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_nom,n,sol.Href(:,1));
% ---- DT propagation
r_transf_orbit_DT  = [SC1.leg1.Xpropad_DT(:,1).*cos(SC1.leg1.Xpropad_DT(:,2)), ...
    SC1.leg1.Xpropad_DT(:,1).*sin(SC1.leg1.Xpropad_DT(:,2)), SC1.leg1.Xpropad_DT(:,3)];
R_cartesian_DT = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_DT,n,sol.Href(:,1));
% ---- AS propagation
r_transf_orbit_AS  = [SC1.leg1.Xpropad_AS(:,1).*cos(SC1.leg1.Xpropad_AS(:,2)), ...
    SC1.leg1.Xpropad_AS(:,1).*sin(SC1.leg1.Xpropad_AS(:,2)), SC1.leg1.Xpropad_AS(:,3)];
R_cartesian_AS = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_AS,n,sol.Href(:,1));

% ----- plots
figure('Name','X comparison')
% plot(SC1.leg1.timead, SC1.leg1.Xad,'Color',colors(1,:),'DisplayName','Xad');
hold on; grid on;
hp1 = plot(SC1.leg1.timead*sim.TU/86400, SC1.leg1.Xpropad_DT,'Color',colors(2,:));
hp2 = plot(SC1.leg1.timead*sim.TU/86400, SC1.leg1.Xpropad_AS,'Color',colors(3,:));
xlabel(' time [d]'); ylabel('dynamics'); % legend([hp1, hp2], 'DT Prop','DT Prop Att');

figure('Name','angles')
plot(SC1.leg1.timead*sim.TU/86400, rad2deg(SC1.leg1.HS.beta),'Color',colors(1,:),'DisplayName','HS el');
hold on; grid on;
plot(SC1.leg1.timead*sim.TU/86400, rad2deg(SC1.leg1.HS.alpha),'Color',colors(1,:),'DisplayName','HS az');
plot(SC1.leg1.timead*sim.TU/86400, rad2deg(SC1.leg1.ANGLES_AS(:,1)),'Color',colors(2,:),'DisplayName','AS az');
plot(SC1.leg1.timead*sim.TU/86400, rad2deg(SC1.leg1.ANGLES_AS(:,2)),'Color',colors(2,:),'DisplayName','HS az');
legend('show'); xlabel(' time [d]'); ylabel('angles [deg]');

figure('Name','Angles rand difference')
plot(SC1.leg1.timead*sim.TU/86400, rad2deg(SC1.leg1.HS.beta - SC1.leg1.ANGLES_AS(:,2)),'Color',colors(1,:),'DisplayName','el diff');
hold on; grid on;
plot(SC1.leg1.timead*sim.TU/86400, rad2deg(SC1.leg1.HS.alpha - SC1.leg1.ANGLES_AS(:,1)),'Color',colors(2,:),'DisplayName','az diff');
legend('show'); xlabel(' time [d]'); ylabel('rand: nominal - AS [deg]');

figure('Name','max err on traj cartesian coordinates')
plot(SC1.leg1.timead*sim.TU/86400, (R_cartesian_DT(:,1) - R_cartesian_AS(:,1))*sim.DU,'Color',colors(1,:),'DisplayName','x diff');
hold on; grid on;
plot(SC1.leg1.timead*sim.TU/86400, (R_cartesian_DT(:,2) - R_cartesian_AS(:,2))*sim.DU,'Color',colors(2,:),'DisplayName','y diff');
plot(SC1.leg1.timead*sim.TU/86400, (R_cartesian_DT(:,3) - R_cartesian_AS(:,3))*sim.DU,'Color',colors(3,:),'DisplayName','z diff');
plot(SC1.leg1.timead*sim.TU/86400, (vecnorm(R_cartesian_DT,2,2) - vecnorm(R_cartesian_AS,2,2))*sim.DU,'Color',colors(4,:),'DisplayName','norm diff');
legend('show'); xlabel(' time [d]'); ylabel('positions [km]');

figure('Name','3D plot orbit')
plot3(R_cartesian_DT(:,1),R_cartesian_DT(:,2),R_cartesian_DT(:,3),'Color',colors(1,:),'DisplayName','DT orbit');
hold on; grid on; axis equal;
plot3(R_cartesian_DT(end,1),R_cartesian_DT(end,2),R_cartesian_DT(end,3),'o','Color',colors(1,:),'DisplayName','DT end');
plot3(R_cartesian_AS(:,1),R_cartesian_AS(:,2),R_cartesian_AS(:,3),'--','Linewidth',3.5,'Color',colors(2,:),'DisplayName','AS orbit');
plot3(R_cartesian_AS(end,1),R_cartesian_AS(end,2),R_cartesian_AS(end,3),'o','Color',colors(2,:),'DisplayName','AS end');
legend('show'); xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');

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
data.fmincon_iter = 500;

SC1.leg2.v_in = v_encounter.astD1;
SC1.leg2.v_end = v_encounter.astA2;
SC1.leg2.r_in = r_encounter.astD1;
SC1.leg2.r_end = r_encounter.astA2;
SC1.leg2.N_rev = sol.Nrev(2);
SC1.leg2.TOF = sol.TOF2_ADIM;
SC1.leg2.Href = sol.Href(:,2);

DT_input = SC1.leg2;
[SC1.leg2.HS, SC1.leg2.ANGLES_AS, SC1.leg2.timead, SC1.leg2.Xad, ...
    SC1.leg2.Xpropad_DT, SC1.leg2.Xpropad_AS] = DT_executable_attitude_sensitivity(DT_input, ...
    sim, data);

%% -- if run with less than 100 points
n = 100;
method = 'linear';
t_new = linspace(SC1.leg2.timead(1),SC1.leg2.timead(end),n)';
SC1.temp.timead       = t_new;
SC1.temp.HS.X = interp1(SC1.leg2.timead,SC1.leg2.HS.X,t_new,method);
SC1.temp.Xpropad_AS = interp1(SC1.leg2.timead,SC1.leg2.Xpropad_AS,t_new,method);
SC1.temp.Xpropad_DT = interp1(SC1.leg2.timead,SC1.leg2.Xpropad_DT,t_new,method);
SC1.temp.HS.beta = interp1(SC1.leg2.timead,SC1.leg2.HS.beta,t_new,method);
SC1.temp.HS.alpha = interp1(SC1.leg2.timead,SC1.leg2.HS.alpha,t_new,method);
SC1.temp.ANGLES_AS = interp1(SC1.leg2.timead,SC1.leg2.ANGLES_AS,t_new,method);

SC1.leg2 = SC1.temp;
SC1 = rmfield(SC1, 'temp');

% ---- cartesian orbit
% ---- nominal propagation
r_transf_orbit_nom  = [SC1.leg2.HS.X(:,1).*cos(SC1.leg2.HS.X(:,2)), ...
    SC1.leg2.HS.X(:,1).*sin(SC1.leg2.HS.X(:,2)), SC1.leg2.HS.X(:,3)];
R_cartesian_nom = rotate_local2ecplitic(r_encounter.astD1,r_transf_orbit_nom,n,sol.Href(:,2));
% ---- DT propagation
r_transf_orbit_DT  = [SC1.leg2.Xpropad_DT(:,1).*cos(SC1.leg2.Xpropad_DT(:,2)), ...
    SC1.leg2.Xpropad_DT(:,1).*sin(SC1.leg2.Xpropad_DT(:,2)), SC1.leg2.Xpropad_DT(:,3)];
R_cartesian_DT = rotate_local2ecplitic(r_encounter.astD1,r_transf_orbit_DT,n,sol.Href(:,2));
% ---- AS propagation
r_transf_orbit_AS  = [SC1.leg2.Xpropad_AS(:,1).*cos(SC1.leg2.Xpropad_AS(:,2)), ...
    SC1.leg2.Xpropad_AS(:,1).*sin(SC1.leg2.Xpropad_AS(:,2)), SC1.leg2.Xpropad_AS(:,3)];
R_cartesian_AS = rotate_local2ecplitic(r_encounter.astD1,r_transf_orbit_AS,n,sol.Href(:,2));

% ----- plots
figure('Name','X comparison')
% plot(SC1.leg1.timead, SC1.leg1.Xad,'Color',colors(1,:),'DisplayName','Xad');
hold on; grid on;
hp1 = plot(SC1.leg2.timead*sim.TU/86400, SC1.leg2.Xpropad_DT,'Color',colors(2,:));
hp2 = plot(SC1.leg2.timead*sim.TU/86400, SC1.leg2.Xpropad_AS,'Color',colors(3,:));
xlabel(' time [d]'); ylabel('dynamics'); % legend([hp1, hp2], 'DT Prop','DT Prop Att');

figure('Name','angles')
plot(SC1.leg2.timead*sim.TU/86400, rad2deg(SC1.leg2.HS.beta),'Color',colors(1,:),'DisplayName','HS el');
hold on; grid on;
plot(SC1.leg2.timead*sim.TU/86400, rad2deg(SC1.leg2.HS.alpha),'Color',colors(1,:),'DisplayName','HS az');
plot(SC1.leg2.timead*sim.TU/86400, rad2deg(SC1.leg2.ANGLES_AS(:,1)),'Color',colors(2,:),'DisplayName','AS az');
plot(SC1.leg2.timead*sim.TU/86400, rad2deg(SC1.leg2.ANGLES_AS(:,2)),'Color',colors(2,:),'DisplayName','HS az');
legend('show'); xlabel(' time [d]'); ylabel('angles [deg]');

figure('Name','Angles rand difference')
plot(SC1.leg2.timead*sim.TU/86400, rad2deg(SC1.leg2.HS.beta - SC1.leg2.ANGLES_AS(:,2)),'Color',colors(1,:),'DisplayName','el diff');
hold on; grid on;
plot(SC1.leg2.timead*sim.TU/86400, rad2deg(SC1.leg2.HS.alpha - SC1.leg2.ANGLES_AS(:,1)),'Color',colors(2,:),'DisplayName','az diff');
legend('show'); xlabel(' time [d]'); ylabel('rand: nominal - AS [deg]');

figure('Name','max err on traj cartesian coordinates')
plot(SC1.leg1.timead*sim.TU/86400, (R_cartesian_DT(:,1) - R_cartesian_AS(:,1))*sim.DU,'Color',colors(1,:),'DisplayName','x diff');
hold on; grid on;
plot(SC1.leg1.timead*sim.TU/86400, (R_cartesian_DT(:,2) - R_cartesian_AS(:,2))*sim.DU,'Color',colors(2,:),'DisplayName','y diff');
plot(SC1.leg1.timead*sim.TU/86400, (R_cartesian_DT(:,3) - R_cartesian_AS(:,3))*sim.DU,'Color',colors(3,:),'DisplayName','z diff');
plot(SC1.leg1.timead*sim.TU/86400, (vecnorm(R_cartesian_DT,2,2) - vecnorm(R_cartesian_AS,2,2))*sim.DU,'Color',colors(4,:),'DisplayName','norm diff');
legend('show'); xlabel(' time [d]'); ylabel('positions [km]');

figure('Name','3D plot orbit')
plot3(R_cartesian_DT(:,1),R_cartesian_DT(:,2),R_cartesian_DT(:,3),'Color',colors(1,:),'DisplayName','DT orbit');
hold on; grid on; axis equal;
plot3(R_cartesian_DT(end,1),R_cartesian_DT(end,2),R_cartesian_DT(end,3),'o','Color',colors(1,:),'DisplayName','DT end');
plot3(R_cartesian_AS(:,1),R_cartesian_AS(:,2),R_cartesian_AS(:,3),'--','Linewidth',3.5,'Color',colors(2,:),'DisplayName','AS orbit');
plot3(R_cartesian_AS(end,1),R_cartesian_AS(end,2),R_cartesian_AS(end,3),'o','Color',colors(2,:),'DisplayName','AS end');
legend('show'); xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
