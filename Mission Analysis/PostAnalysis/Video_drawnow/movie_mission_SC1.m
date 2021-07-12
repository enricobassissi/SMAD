%% MOVIE PLOT OF THE WHOLE MISSION SC1
%{
and the full orbit of inner solar system in that same 
time of the mission
%}
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
str_path=split(pwd, 'PostAnalysis\Video_drawnow');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
dt_path=string(str_path(1))+'TrajOptimisation\direct transcript\follower';
addpath(genpath(dt_path));
fun_path=string(str_path(1))+'TrajOptimisation\direct transcript\functions';
addpath(genpath(fun_path));

%% load and textures
load('Mission_160kg_dry.mat')
load('MA5_160kg_dry.mat')
load('data_elements_matrix_44_63_2SC.mat')
% load('new_miss_vid_data3.mat')

%% --------------------------------------------------------------------- %%
%% --------------------------------- SC1 ------------------------------- %%
%% --------------------------------------------------------------------- %%
%% the time vector
% MA.SC1.uniform.time = [SC1.mjd2000_dep+SC1.leg1.timead.*sim.TU/86400; SC1.coasting.leg1.time;
%     SC1.coasting.leg1.time(end)+SC1.leg2.timead.*sim.TU/86400; SC1.coasting.leg2.time;]; % mjd 2000
time_vect = MA.SC1.mjd2000_dep+MA.SC1.uniform.time/86400;
err_1st_tof = time_vect(101)-time_vect(100)
new_tof1 = time_vect(100)-time_vect(1)
new_ct1 = SC1.CT1+time_vect(101)-time_vect(100)
err_2nd_tof = time_vect(301)-time_vect(300)
new_tof2 = time_vect(300)-time_vect(201)

% time_vect2 = MA.SC2.mjd2000_dep+MA.SC2.uniform.time/86400;
% time_vect2(101)-time_vect2(100)
% time_vect2(301)-time_vect2(300)
% -- length
n = length(time_vect);

%% sc trajectory
% --- new encounters
MJD01_dim = SC1.mjd2000_dep;
[kep_EA,~] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,sim.mu_dim);
r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;
new_r_encounter.EA = r_EA;
new_v_encounter.EA = v_EA;

MJDA1_dim = (SC1.mjd2000_dep+new_tof1);
[kep_ast_A1] = uNEO3(MJDA1_dim,SC1.asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[rA1, vA1] = sv_from_coe(kep_ast_A1,sim.mu_dim); % km, km/s
rA1 = rA1/sim.DU;
vA1 = vA1/sim.DU*sim.TU;
new_r_encounter.astA1 = rA1;
new_v_encounter.astA1 = vA1;

MJDD1_dim = (SC1.mjd2000_dep+new_tof1+new_ct1);
[kep_ast_D1] = uNEO3(MJDD1_dim,SC1.asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[rD1, vD1] = sv_from_coe(kep_ast_D1,sim.mu_dim); % km, km/s
rD1 = rD1/sim.DU;
vD1 = vD1/sim.DU*sim.TU;
new_r_encounter.astD1 = rD1;
new_v_encounter.astD1 = vD1;

MJDA2_dim = (SC1.mjd2000_dep+new_tof1+new_ct1+new_tof2);
[kep_ast_A2] = uNEO3(MJDA2_dim,SC1.asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[rA2, vA2] = sv_from_coe(kep_ast_A2,sim.mu_dim); % km, km/s
rA2 = rA2/sim.DU;
vA2 = vA2/sim.DU*sim.TU;
new_r_encounter.astA2 = rA2;
new_v_encounter.astA2 = vA2;

MJDD2_dim = (SC1.mjd2000_dep+new_tof1+new_ct1+new_tof2+new_ct1);
[kep_ast_D2] = uNEO3(MJDD2_dim,SC1.asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[rD2, vD2] = sv_from_coe(kep_ast_D2,sim.mu_dim); % km, km/s
rD2 = rD2/sim.DU;
vD2 = vD2/sim.DU*sim.TU;
new_r_encounter.astD2 = rD2;
new_v_encounter.astD2 = vD2;

%% DT - SC 1 - LEG 1
data.Tmax = 0.02; % max thrust available for EGI
data.Pmax = 550; % max power related to 20 mN thrust on qitetiq t5
data.Isp = 2900/sim.TU; % for 20 mN the Isp of Qinetiq T5 is 2900s
data.n_int = 20; % discretisation selected
data.angle_inplane_panels_max = 20; % degrees
data.MU = sol.M_start_SC1_leg1; % mass adimensionalisation on wet mass
SC1.asteroid_1 = sol.asteroid_1;
data.muS = sim.mu_dim;   
data.ThrustModulationFlag = 1;

v_launcher = sol.v_inf_magn/sim.DU*sim.TU*[cos(sol.el)*cos(sol.az); cos(sol.el)*sin(sol.az); sin(sol.el)];
v_dep = new_v_encounter.EA + v_launcher;  %if parabolic escape (v_extra = 0)
SC1.leg1.v_in = v_dep;
SC1.leg1.v_end = new_v_encounter.astA1;
SC1.leg1.r_in = new_r_encounter.EA;
SC1.leg1.r_end = new_r_encounter.astA1;
SC1.leg1.N_rev = sol.Nrev(1);
SC1.leg1.TOF = new_tof1*86400/sim.TU;
SC1.leg1.Href = sol.Href(:,1);

DT_input = SC1.leg1;
tic
[SC1.leg1.HS, SC1.leg1.RES, SC1.leg1.timead, SC1.leg1.Xad, SC1.leg1.Xpropad, ...
    SC1.leg1.Xpropad_DT] = DT_executable(DT_input, sim, data);
el_time.DT_SC1_leg1 = toc/60; % minutes

SC1.leg1.PROPULSION = nikita(SC1.leg1, data, sim, colors);
SC1.leg1.EPS = marco_elia(SC1.leg1,new_r_encounter.EA,data,sim,SC1.leg1.Href,colors);

%% -- if run with less than 100 points
data.n_int = 100;
[SC1.leg1] = interpoliamo_porc(SC1.leg1,sim.n_sol);
SC1.leg1.v_in = v_dep;
SC1.leg1.v_end = new_v_encounter.astA1;
SC1.leg1.r_in = new_r_encounter.EA;
SC1.leg1.r_end = new_r_encounter.astA1;
SC1.leg1.N_rev = sol.Nrev(1);
SC1.leg1.TOF = new_tof1*86400/sim.TU;
SC1.leg1.Href = sol.Href(:,1);
SC1.mjd2000_dep = sol.departure_mjd2000;

SC1.leg1.PROPULSION = nikita(SC1.leg1, data, sim, colors);
SC1.leg1.EPS = marco_elia(SC1.leg1,new_r_encounter.EA,data,sim,SC1.leg1.Href,colors);

% -- new mass after DT calculations
SC1.leg1.mass_start = SC1.leg1.HS.X(1,7);
SC1.leg1.mass_end = SC1.leg1.HS.X(end,7);
SC1.leg1.mass_depleted = SC1.leg1.mass_start - SC1.leg1.mass_end;
SC1.leg1.mass_fraction = SC1.leg1.mass_depleted/SC1.leg1.mass_start;

% -- coasting between the two legs
[SC1.coasting.leg1] = relative_position_coasting_stuff(SC1.asteroid_1,SC1.mjd2000_dep,...
    new_tof1,new_ct1,sim.n_sol,EphData,sim,colors);%leg1.timead(end)*sim.TU/86400

%% DT - SC 1 - LEG 2
data.Tmax = 0.015; % max thrust available for EGI
data.Pmax = 450; % max power related to 15 mN thrust on qitetiq t5
data.Isp = 2750/sim.TU; % for 15 mN the Isp of Qinetiq T5 is 2750s
data.ThrustModulationFlag = 2;
data.angle_inplane_panels_max = 10; % degrees
data.n_int = 50;
% data.MU = sol.M_start_SC1_leg2; % mass adimensionalisation on wet mass
data.MU = SC1.leg1.mass_end - sim.M_pods; % mass adimensionalisation on wet mass
SC1.asteroid_2 = sol.asteroid_2;

SC1.leg2.v_in = new_v_encounter.astD1;
SC1.leg2.v_end = new_v_encounter.astA2;
SC1.leg2.r_in = new_r_encounter.astD1;
SC1.leg2.r_end = new_r_encounter.astA2;
SC1.leg2.N_rev = sol.Nrev(2);
SC1.leg2.TOF = new_tof2*86400/sim.TU;
SC1.leg2.Href = sol.Href(:,2);

DT_input = SC1.leg2;
tic
[SC1.leg2.HS, SC1.leg2.RES, SC1.leg2.timead, SC1.leg2.Xad, SC1.leg2.Xpropad, ...
    SC1.leg2.Xpropad_DT] = DT_executable(DT_input, sim, data);
el_time.DT_SC1_leg2 = toc/60; % minutes

SC1.leg2.PROPULSION = nikita(SC1.leg2, data, sim, colors);
SC1.leg2.EPS = marco_elia(SC1.leg2,new_r_encounter.astD1,data,sim,SC1.leg2.Href,colors);

%% -- if run with less than 100 points
data.n_int = 100;
[SC1.leg2] = interpoliamo_porc(SC1.leg2,data.n_int);
SC1.leg2.v_in = new_v_encounter.astD1;
SC1.leg2.v_end = new_v_encounter.astA2;
SC1.leg2.r_in = new_r_encounter.astD1;
SC1.leg2.r_end = new_r_encounter.astA2;
SC1.leg2.N_rev = sol.Nrev(2);
SC1.leg2.TOF = new_tof2*86400/sim.TU;
SC1.leg2.Href = sol.Href(:,2);
SC1.mjd2000_dep = sol.departure_mjd2000;

SC1.leg2.PROPULSION = nikita(SC1.leg2, data, sim, colors);
SC1.leg2.EPS = marco_elia(SC1.leg2,r_encounter.astD1,data,sim,SC1.leg2.Href,colors);

% -- new mass after DT calculations
SC1.leg2.mass_start = SC1.leg2.HS.X(1,7);
SC1.leg2.mass_end = SC1.leg2.HS.X(end,7);
SC1.leg2.mass_depleted = SC1.leg2.mass_start - SC1.leg2.mass_end;
SC1.leg2.mass_fraction = SC1.leg2.mass_depleted/SC1.leg2.mass_start;

% -- coasting between the two legs
SC1.mjd2000_dep_ast1 = SC1.mjd2000_dep+new_tof1+new_ct1;
[SC1.coasting.leg2] = relative_position_coasting_stuff(SC1.asteroid_2,SC1.mjd2000_dep_ast1,...
    new_tof2,new_ct1,sim.n_sol,EphData,sim,colors); %SC1.leg2.timead(end)*sim.TU/86400

%% -- building of a uniform position, time, thrust vector for SC1
SC1.uniform.R = [SC1.leg1.EPS.R_cartesian; SC1.coasting.leg1.r_ast; 
    SC1.leg2.EPS.R_cartesian; SC1.coasting.leg2.r_ast;]./sim.DU; % AU
figure()
% plot3(SC1.uniform.R(:,1),SC1.uniform.R(:,2),SC1.uniform.R(:,3));
plot3(SC1.uniform.R(1:sim.n_sol,1),SC1.uniform.R(1:sim.n_sol,2),SC1.uniform.R(1:sim.n_sol,3),...
    'Color',colors(1,:),'DisplayName','Cruise - Leg 1');
axis equal; grid on; hold on;
plot3(SC1.uniform.R(sim.n_sol+1:2*sim.n_sol,1),SC1.uniform.R(sim.n_sol+1:2*sim.n_sol,2),...
    SC1.uniform.R(sim.n_sol+1:2*sim.n_sol,3),...
    'Color',colors(2,:),'DisplayName','Rendezvous - 2020VV');
plot3(SC1.uniform.R(2*sim.n_sol+1:3*sim.n_sol,1),SC1.uniform.R(2*sim.n_sol+1:3*sim.n_sol,2),...
    SC1.uniform.R(2*sim.n_sol+1:3*sim.n_sol,3),...
    'Color',colors(3,:),'DisplayName','Cruise - Leg 2');
plot3(SC1.uniform.R(3*sim.n_sol+1:4*sim.n_sol,1),SC1.uniform.R(3*sim.n_sol+1:4*sim.n_sol,2),...
    SC1.uniform.R(3*sim.n_sol+1:4*sim.n_sol,3),...
    'Color',colors(4,:),'DisplayName','Rendezvous - 2009TD17');
plot3(new_r_encounter.EA(1), new_r_encounter.EA(2), new_r_encounter.EA(3),'*',...
    'Color',colors(8,:),'DisplayName','Earth Dep');
plot3(new_r_encounter.astA1(1), new_r_encounter.astA1(2), new_r_encounter.astA1(3),'^',...
    'DisplayName','2020VV Arr');
plot3(new_r_encounter.astD1(1), new_r_encounter.astD1(2), new_r_encounter.astD1(3),'*',...
    'DisplayName','2020VV Dep');
plot3(new_r_encounter.astA2(1), new_r_encounter.astA2(2), new_r_encounter.astA2(3),'^',...
    'DisplayName','2009TD17 Arr');
plot3(new_r_encounter.astD2(1), new_r_encounter.astD2(2), new_r_encounter.astD2(3),'*',...
    'DisplayName','End Mission');
hsun = plot3(0,0,0,'*','Color',colors(4,:));
hsun.Annotation.LegendInformation.IconDisplayStyle = 'off';

[~, h_earth_whole] = plot_object_orbit(new_time_vect(1),'earth',365,sim,data,100,colors,8);
[~, h_mars_whole] = plot_object_orbit(new_time_vect(1),'mars',687,sim,data,100,colors,6);
h_earth_whole.Annotation.LegendInformation.IconDisplayStyle = 'off';
h_mars_whole.Annotation.LegendInformation.IconDisplayStyle = 'off';
h_earth_whole.Color(4) = 0.2;
h_mars_whole.Color(4) = 0.2;
legend('show','Location','northeastoutside');
view(2)
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');

SC1.uniform.time = [SC1.mjd2000_dep+SC1.leg1.timead.*sim.TU/86400; SC1.coasting.leg1.time;
    SC1.coasting.leg1.time(end)+SC1.leg2.timead.*sim.TU/86400; SC1.coasting.leg2.time;]; % mjd 2000

SC1.uniform.Thrust = [SC1.leg1.EPS.T_transf_orbit; zeros(sim.n_sol,3);
    SC1.leg2.EPS.T_transf_orbit; zeros(sim.n_sol,3)]; % N
figure()
plot(SC1.uniform.time, SC1.uniform.Thrust)
xlabel('t [mjd2000]'); ylabel('T [N]'); legend('Tx','Ty','Tz');
figure()
plot(SC1.uniform.time, vecnorm(SC1.uniform.Thrust,2,2))
xlabel('t [mjd2000]'); ylabel('T magnitude [N]');

%% time interpolation
time_vect = SC1.uniform.time;
time_vect(100) = (time_vect(101)+time_vect(99))/2;
time_vect(201) = (time_vect(202)+time_vect(200))/2;
time_vect(299) = time_vect(298)+(time_vect(301)-time_vect(298))/3;
time_vect(300) = time_vect(298)+(time_vect(301)-time_vect(298))*2/3;

dt_cruise = time_vect(2)-time_vect(1);
dt_rv = time_vect(103)-time_vect(102);
number_of_new_time_intervals = round((time_vect(end)-time_vect(1))/dt_rv);
new_time_vect = linspace(time_vect(1),time_vect(end),number_of_new_time_intervals)';

n = length(new_time_vect);

%% full orbit propagations
% --- same time of propagations between planets/asteroids and sc
for k=1:n
    % -- mercury
% 	[kepME,ksun] = uplanet(time_vect(k), 1);
% 	[rME(k,:), vME(k,:)] = sv_from_coe(kepME,ksun);  
% 	V1(k,:)=vME;
    % -- venus
% 	[kepVE,ksun] = uplanet(time_vect(k), 2);
% 	[rVE(k,:), vVE] = sv_from_coe(kepVE,ksun);  
% 	V2(k,:)=vVE;
    % -- Earth
	[kepE,ksun] = uplanet(new_time_vect(k), 3);
	[rEA(k,:), vEA] = sv_from_coe(kepE,ksun);  
% 	V3(k,:)=vE;
    % -- mars
	[kepMA,ksun] = uplanet(new_time_vect(k), 4);
	[rMA(k,:), vMA] = sv_from_coe(kepMA,ksun);  
% 	V4(k,:)=vMA;
    
    % --- ast 1
    [kep_ast1] = uNEO3(new_time_vect(k),MA.SC1.asteroid_1,data);
    [r_ast1(k,1:3),v_ast1] = sv_from_coe(kep_ast1,ksun); % km, km/s
    % --- ast 2
    [kep_ast2] = uNEO3(new_time_vect(k),MA.SC1.asteroid_2,data);
    [r_ast2(k,1:3),v_ast2] = sv_from_coe(kep_ast2,ksun); % km, km/s
end
clearvars k

% -- adimensionalise
% rME=rME/sim.DU; rVE=rVE/sim.DU; 
rEA=rEA/sim.DU; rMA=rMA/sim.DU;
r_ast1 = r_ast1/sim.DU; r_ast2 = r_ast2/sim.DU;

%% sc trajectory interpolation
[az_SC1_sph_pre,elev_SC1_sph_pre,r_SC1_sph_pre] = cart2sph(SC1.uniform.R(:,1),SC1.uniform.R(:,2),SC1.uniform.R(:,3));
az_SC1_sph_post = smooth(interp1(time_vect,smooth(unwrap(az_SC1_sph_pre),6),new_time_vect),20);
elev_SC1_sph_post = interp1(time_vect,smooth(elev_SC1_sph_pre,5),new_time_vect);
r_SC1_sph_post = interp1(time_vect,smooth(r_SC1_sph_pre,5),new_time_vect);
[xpl,ypl,zpl] = sph2cart(wrapToPi(az_SC1_sph_post),elev_SC1_sph_post,r_SC1_sph_post);
% yplot = 0;
% yplot = SC1.uniform.R;
yplot = [xpl,ypl,zpl];

figure()
hold on
plot(time_vect,unwrap(az_SC1_sph_pre))
plot(new_time_vect,az_SC1_sph_post)
% plot(time_vect,smooth(unwrap(az_SC1_sph_pre),10))
% plot(time_vect,smooth(unwrap(az_SC1_sph_pre),6))
title('az')

figure()
plot(time_vect,elev_SC1_sph_pre)
hold on
plot(new_time_vect,elev_SC1_sph_post)
title('elev')

figure()
plot(time_vect,r_SC1_sph_pre)
hold on
plot(new_time_vect,r_SC1_sph_post)
title('r')

figure()
plot(time_vect,SC1.uniform.R(:,1))
hold on
plot(new_time_vect,xpl)
title('x')

figure()
plot(time_vect,SC1.uniform.R(:,2))
hold on
plot(new_time_vect,ypl)
title('y')

figure()
plot(time_vect,SC1.uniform.R(:,3))
hold on
plot(new_time_vect,zpl)
title('z')

%% Thrust stuff
T_SC1 = [SC1.leg1.HS.T;zeros(sim.n_sol,1);SC1.leg2.HS.T;zeros(sim.n_sol,1)];
alpha_SC1 = [SC1.leg1.HS.alpha;zeros(sim.n_sol,1);SC1.leg2.HS.alpha;zeros(sim.n_sol,1)];
beta_SC1 = [SC1.leg1.HS.beta;zeros(sim.n_sol,1);SC1.leg2.HS.beta;zeros(sim.n_sol,1)];

figure()
plot(time_vect,T_SC1)
xlabel('t [mjd2000]'); ylabel('T [N]'); 
figure()
plot(time_vect,alpha_SC1)
xlabel('t [mjd2000]'); ylabel('T [N]'); 
figure()
plot(time_vect,beta_SC1)
xlabel('t [mjd2000]'); ylabel('T [N]'); 

% --- reinterpolate
T_SC1_post = interp1(time_vect,T_SC1,new_time_vect);
alpha_SC1_post = interp1(time_vect,alpha_SC1,new_time_vect);
beta_SC1_post = interp1(time_vect,beta_SC1,new_time_vect);

figure()
plot(time_vect,T_SC1)
hold on
plot(new_time_vect,T_SC1_post)
xlabel('t [mjd2000]'); ylabel('T [N]');

figure()
plot(time_vect,beta_SC1)
hold on
plot(new_time_vect,beta_SC1_post)
xlabel('t [mjd2000]'); ylabel('beta[rad]'); 
figure()
plot(time_vect,alpha_SC1)
hold on
plot(new_time_vect,alpha_SC1_post)
xlabel('t [mjd2000]'); ylabel('alpha[rad]'); 

%% video stuff SC1
load('data_before_video_sc1.mat')

% --- limits of the camera pov
lim_max_helio = max(vecnorm(rMA,2,2));
% rel_sc_ast1 = yplot-r_ast1;
% % lim_max_ast1 = max(vecnorm(rel_sc_ast1,2,2));
% lim_max_ast1 =0.15;
% min_rel_sc_ast1_non_zero = 10;
% for i = 1:n
%     if norm(rel_sc_ast1(i,:)) > 0
%         temp_nz = norm(rel_sc_ast1(i,:));
%         if temp_nz < min_rel_sc_ast1_non_zero
%             min_rel_sc_ast1_non_zero = temp_nz;
%         end
%     end
% end
% clearvars i
% rel_sc_ast2 = yplot-r_ast2;
% min_rel_sc_ast2_non_zero = 10;
% lim_max_ast2 = max(vecnorm(rel_sc_ast2,2,2));
% for i = 1:n
%     if norm(rel_sc_ast2(i,:)) > 0
%         temp_nz = norm(rel_sc_ast2(i,:));
%         if temp_nz < min_rel_sc_ast2_non_zero
%             min_rel_sc_ast2_non_zero = temp_nz;
%         end
%     end
% end
% clearvars i

tic
% --- movie
figure('Name','Mission Video SC1'); %,'WindowState','maximized'
set(gca,'nextplot','replacechildren');
set(gcf, 'Position', get(0, 'Screensize'));
v = VideoWriter('mission_SC1_1','MPEG-4');
% v = VideoWriter('mission4', 'Uncompressed AVI');
v.Quality = 100;
open(v);  

% set(gca, 'CameraPosition', [0 lim_max lim_max]);
for i=1:2:n % i=1:6:n % i=1:389:n
    % ---- start orbit subplot
    subplot(4,3,[1,2,4,5,7,8,10,11])
    % --- whole orbits
    axis equal; grid on; hold on;
    xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
    [~, h_ast1_whole] = plot_object_orbit(new_time_vect(1),'asteroid',1.5*365,sim,data,100,colors,3,MA.SC1.asteroid_1);
    [~, h_ast2_whole] = plot_object_orbit(new_time_vect(1),'asteroid',1.5*365,sim,data,100,colors,4,MA.SC1.asteroid_2);
    [~, h_earth_whole] = plot_object_orbit(new_time_vect(1),'earth',365,sim,data,100,colors,8);
    [~, h_mars_whole] = plot_object_orbit(new_time_vect(1),'mars',687,sim,data,100,colors,6);
    h_sun = plot3(0,0,0,'*','Color',colors(4,:),'MarkerSize',5); % sun
    h_ast1_whole.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_ast2_whole.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_earth_whole.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_mars_whole.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_sun.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_ast1_whole.Color(4) = 0.3;
    h_ast2_whole.Color(4) = 0.3;
    h_earth_whole.Color(4) = 0.2;
    h_mars_whole.Color(4) = 0.2;
    
    % --- highlight the actual planet position
    [~, h_ast1_highlight] = plot_object_orbit(new_time_vect(i),'asteroid',0.25*365,sim,data,30,colors,3,MA.SC1.asteroid_1);
    [~, h_ast2_highlight] = plot_object_orbit(new_time_vect(i),'asteroid',0.25*365,sim,data,30,colors,4,MA.SC1.asteroid_2);
    [~, h_earth_highlight] = plot_object_orbit(new_time_vect(i),'earth',0.2*365,sim,data,30,colors,8);
    [~, h_mars_highlight] = plot_object_orbit(new_time_vect(i),'mars',0.3*365,sim,data,30,colors,6);
    h_ast1_highlight.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_ast2_highlight.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_earth_highlight.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_mars_highlight.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_ast1_highlight.Color(4) = 0.5;
    h_ast2_highlight.Color(4) = 0.5;
    h_earth_highlight.Color(4) = 0.3;
    h_mars_highlight.Color(4) = 0.3;
    
    [~, h_ast1_highlight2] = plot_object_orbit(new_time_vect(i),'asteroid',0.08*365,sim,data,10,colors,3,MA.SC1.asteroid_1);
    [~, h_ast2_highlight2] = plot_object_orbit(new_time_vect(i),'asteroid',0.08*365,sim,data,10,colors,4,MA.SC1.asteroid_2);
    [~, h_earth_highlight2] = plot_object_orbit(new_time_vect(i),'earth',0.06*365,sim,data,10,colors,8);
    [~, h_mars_highlight2] = plot_object_orbit(new_time_vect(i),'mars',0.1*365,sim,data,10,colors,6);
    h_ast1_highlight2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_ast2_highlight2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_earth_highlight2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_mars_highlight2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_ast1_highlight2.Color(4) = 1;
    h_ast2_highlight2.Color(4) = 1;
    h_earth_highlight2.Color(4) = 0.6;
    h_mars_highlight2.Color(4) = 0.6;
    
    % --- the coma of the spacecraft trajectory
    if i > 41
        k = 40;
        h_coma_sc = plot3( yplot(i-k:i,1), yplot(i-k:i,2), yplot(i-k:i,3),...
                '--','Color',colors(1,:),'LineWidth',1);
        h_coma_sc.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h_coma_sc.Color(4) = 0.5;
    end
    if i > 21
        k=20;
        h_coma_sc2 = plot3( yplot(i-k:i,1), yplot(i-k:i,2), yplot(i-k:i,3),...
                '--','Color',colors(1,:),'LineWidth',1);
        h_coma_sc2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h_coma_sc2.Color(4) = 1;
    end

    % --- normal condition, nominal trajectory heliocentric, point by point
    h0 = plot3( yplot(i,1), yplot(i,2), yplot(i,3),'s','Color',colors(1,:),'MarkerSize',9,...
        'DisplayName','SC1');
    h3 = plot3( rEA(i,1), rEA(i,2), rEA(i,3),'o','Color',colors(8,:),'MarkerSize',4,...
        'LineWidth',1.5,'DisplayName','Earth');
    h4 = plot3( rMA(i,1), rMA(i,2), rMA(i,3),'o','Color',colors(6,:),'MarkerSize',4,...
        'LineWidth',1.5,'DisplayName','Mars');
    hast1 = plot3( r_ast1(i,1), r_ast1(i,2), r_ast1(i,3),'*','Color',colors(3,:),...
        'MarkerSize',5,'DisplayName',MA.SC1.asteroid_1);
    hast2 = plot3( r_ast2(i,1), r_ast2(i,2), r_ast2(i,3),'*','Color',colors(4,:),...
        'MarkerSize',5,'DisplayName',MA.SC1.asteroid_2);

    legend('show','Location','northwest')
    
    date_view=mjd20002date(new_time_vect(i));
    title (sprintf('Mission Date %d %d %d', date_view(1:3)));
%     set(get(gca,'title'),'Position',[15 13 1.00011])

    view([lim_max_helio lim_max_helio lim_max_helio])  
    xlim([-lim_max_helio lim_max_helio])
    ylim([-lim_max_helio lim_max_helio])
    zlim([-lim_max_helio lim_max_helio])

    % -- textbox inside the orbit plots
    if new_time_vect(i) < MJDA1_dim
        cs1 = '\begin{tabular}{l}  \textbf{Phase} \\ Cruise - Leg 1 \\ Rendezvous - 2020VV \\ Cruise - Leg 2 \\ Rendezvous - 2009TD17 \end{tabular}';
        cs2 = {'\textbf{Mission Time [d]}'; ['T+', num2str(round(new_time_vect(i)-SC1.mjd2000_dep))] ;
            ['T-', num2str(round(MJDA1_dim-new_time_vect(i)))];
            '-'; '-'};
        h_textbox1 = annotation('textbox', [0.28, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        h_textbox2 = annotation('textbox', [0.40, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        set(h_textbox1, 'String', cs1);
        set(h_textbox2, 'String', cs2);
    elseif new_time_vect(i) >= MJDA1_dim && new_time_vect(i) < MJDD1_dim
        cs1 = '\begin{tabular}{l}  \textbf{Phase} \\ Cruise - Leg 1 \\ Rendezvous - 2020VV \\ Cruise - Leg 2 \\ Rendezvous - 2009TD17 \end{tabular}';
        cs2 = {'\textbf{Mission Time [d]}'; 'Completed' ; 
            ['T+', num2str(round(new_time_vect(i)-MJDA1_dim))];
            ['T-', num2str(round(MJDD1_dim-new_time_vect(i)))]; '-'};
        h_textbox1 = annotation('textbox', [0.28, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        h_textbox2 = annotation('textbox', [0.40, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        set(h_textbox1, 'String', cs1);
        set(h_textbox2, 'String', cs2);
    elseif new_time_vect(i) >= MJDD1_dim && new_time_vect(i) < MJDA2_dim
        cs1 = '\begin{tabular}{l}  \textbf{Phase} \\ Cruise - Leg 1 \\ Rendezvous - 2020VV \\ Cruise - Leg 2 \\ Rendezvous - 2009TD17 \end{tabular}';
        cs2 = {'\textbf{Mission Time [d]}'; 'Completed' ; 'Completed';
            ['T+', num2str(round(new_time_vect(i)-MJDD1_dim))]; 
            ['T-', num2str(round(MJDA2_dim-new_time_vect(i)))]};
        h_textbox1 = annotation('textbox', [0.28, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        h_textbox2 = annotation('textbox', [0.40, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        set(h_textbox1, 'String', cs1);
        set(h_textbox2, 'String', cs2);
	elseif new_time_vect(i) >= MJDA2_dim
        cs1 = '\begin{tabular}{l}  \textbf{Phase} \\ Cruise - Leg 1 \\ Rendezvous - 2020VV \\ Cruise - Leg 2 \\ Rendezvous - 2009TD17 \end{tabular}';
        cs2 = {'\textbf{Mission Time [d]}'; 'Completed' ; 'Completed'; 'Completed';
             ['T+', num2str(round(new_time_vect(i)-MJDA2_dim))]};
        h_textbox1 = annotation('textbox', [0.28, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        h_textbox2 = annotation('textbox', [0.40, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        set(h_textbox1, 'String', cs1);
        set(h_textbox2, 'String', cs2);
    end
    
    drawnow

    % ------ thrust profile subplot
    subplot(4,3,[3,6])
    h_plot_subplot3 = plot(new_time_vect-SC1.mjd2000_dep,T_SC1_post,'Color',colors(1,:));
    hold on; grid on;
    h_plot_subplot3_point = plot(new_time_vect(i)-SC1.mjd2000_dep,T_SC1_post(i),'o','Color',colors(2,:));
%     xlabel('T+dep [d]'); 
    ylabel('T [N]'); 
    xlim([min(new_time_vect-SC1.mjd2000_dep) max(new_time_vect-SC1.mjd2000_dep)])
    xticks([0 500 1000 1500 2000])
    ylim([min(T_SC1_post) max(T_SC1_post)])
    
    h_rvd1_region = fill([MJDA1_dim-SC1.mjd2000_dep,MJDD1_dim-SC1.mjd2000_dep,...
        MJDD1_dim-SC1.mjd2000_dep,MJDA1_dim-SC1.mjd2000_dep],...
        [min(T_SC1_post),min(T_SC1_post),...
        max(T_SC1_post),max(T_SC1_post)],'red');
    h_rvd1_region.FaceAlpha=0.3;
    h_rvd1_region.FaceColor=colors(10,:);
    h_rvd1_region.EdgeColor='none';
    
    h_rvd2_region = fill([MJDA2_dim-SC1.mjd2000_dep,MJDD2_dim-SC1.mjd2000_dep,...
        MJDD2_dim-SC1.mjd2000_dep,MJDA2_dim-SC1.mjd2000_dep],...
        [min(T_SC1_post),min(T_SC1_post),...
        max(T_SC1_post),max(T_SC1_post)],'red');
    h_rvd2_region.FaceAlpha=0.3;
    h_rvd2_region.FaceColor=colors(10,:);
    h_rvd2_region.EdgeColor='none';
    
    legend([h_rvd1_region,h_plot_subplot3,h_plot_subplot3_point],...
        {'Rendezvous Periods','Thrust Profile','Mission Moment'},...
        'Location','northoutside');
    h_rvd2_region.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     h_plot_subplot3.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     h_plot_subplot3_point.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    drawnow

    subplot(4,3,9)
    h_plot_subplot61 = plot(new_time_vect-SC1.mjd2000_dep,wrapTo180(rad2deg(alpha_SC1_post)),'Color',colors(1,:));
    hold on; grid on;
    h_plot_subplot61_point = plot(new_time_vect(i)-SC1.mjd2000_dep,wrapTo180(rad2deg(alpha_SC1_post(i))),'o','Color',colors(2,:));
%     xlabel('T+dep [d]'); 
    ylabel('Azimuth [deg]'); 
    xlim([min(new_time_vect-SC1.mjd2000_dep) max(new_time_vect-SC1.mjd2000_dep)])
    xticks([0 500 1000 1500 2000])
    ylim([min(rad2deg(alpha_SC1_post)) max(rad2deg(alpha_SC1_post))])
    
    h_rvd3_region = fill([MJDA1_dim-SC1.mjd2000_dep,MJDD1_dim-SC1.mjd2000_dep,...
        MJDD1_dim-SC1.mjd2000_dep,MJDA1_dim-SC1.mjd2000_dep],...
        [min(rad2deg(alpha_SC1_post)),min(rad2deg(alpha_SC1_post)),...
        max(rad2deg(alpha_SC1_post)),max(rad2deg(alpha_SC1_post))],'red');
    h_rvd3_region.FaceAlpha=0.3;
    h_rvd3_region.FaceColor=colors(10,:);
    h_rvd3_region.EdgeColor='none';
    
    h_rvd4_region = fill([MJDA2_dim-SC1.mjd2000_dep,MJDD2_dim-SC1.mjd2000_dep,...
        MJDD2_dim-SC1.mjd2000_dep,MJDA2_dim-SC1.mjd2000_dep],...
        [min(rad2deg(alpha_SC1_post)),min(rad2deg(alpha_SC1_post)),...
        max(rad2deg(alpha_SC1_post)),max(rad2deg(alpha_SC1_post))],'red');
    h_rvd4_region.FaceAlpha=0.3;
    h_rvd4_region.FaceColor=colors(10,:);
    h_rvd4_region.EdgeColor='none';
    
    h_plot_subplot61.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_plot_subplot61_point.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    subplot(4,3,12)
    h_plot_subplot62 = plot(new_time_vect-SC1.mjd2000_dep,wrapTo180(rad2deg(beta_SC1_post)),'Color',colors(1,:));
    hold on; grid on;
    h_plot_subplot62_point = plot(new_time_vect(i)-SC1.mjd2000_dep,wrapTo180(rad2deg(beta_SC1_post(i))),'o','Color',colors(2,:));
    xlabel('T+dep [d]'); ylabel('Elevation [deg]'); 
    xlim([min(new_time_vect-SC1.mjd2000_dep) max(new_time_vect-SC1.mjd2000_dep)])
    xticks([0 500 1000 1500 2000])
    ylim([min(rad2deg(beta_SC1_post)) max(rad2deg(beta_SC1_post))])
    
    h_rvd5_region = fill([MJDA1_dim-SC1.mjd2000_dep,MJDD1_dim-SC1.mjd2000_dep,...
        MJDD1_dim-SC1.mjd2000_dep,MJDA1_dim-SC1.mjd2000_dep],...
        [min(rad2deg(beta_SC1_post)),min(rad2deg(beta_SC1_post)),...
        max(rad2deg(beta_SC1_post)),max(rad2deg(beta_SC1_post))],'red');
    h_rvd5_region.FaceAlpha=0.3;
    h_rvd5_region.FaceColor=colors(10,:);
    h_rvd5_region.EdgeColor='none';
    
    h_rvd6_region = fill([MJDA2_dim-SC1.mjd2000_dep,MJDD2_dim-SC1.mjd2000_dep,...
        MJDD2_dim-SC1.mjd2000_dep,MJDA2_dim-SC1.mjd2000_dep],...
        [min(rad2deg(beta_SC1_post)),min(rad2deg(beta_SC1_post)),...
        max(rad2deg(beta_SC1_post)),max(rad2deg(beta_SC1_post))],'red');
    h_rvd6_region.FaceAlpha=0.3;
    h_rvd6_region.FaceColor=colors(10,:);
    h_rvd6_region.EdgeColor='none';
    
    h_plot_subplot62.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_plot_subplot62_point.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    drawnow
    
    % --- get the frame in the video
    frame = getframe(gcf);
    writeVideo(v,frame);
    
    % --- delete stuff
    % -- orbits
    delete(h0);
    delete(h3);
    delete(h4);
    delete(hast1);
    delete(hast2);
    % -- highlight
    delete(h_ast1_highlight);
    delete(h_ast2_highlight);
    delete(h_earth_highlight);
    delete(h_mars_highlight);
    delete(h_ast1_highlight2);
    delete(h_ast2_highlight2);
    delete(h_earth_highlight2);
    delete(h_mars_highlight2);
    if i>41
        delete(h_coma_sc);
    end
    if i>21
        delete(h_coma_sc2);
    end
    % -- thrust
    delete(h_plot_subplot3_point);
    % -- angles
    delete(h_plot_subplot61_point);
    delete(h_plot_subplot62_point);
    % -- text box
    delete(h_textbox1);
    delete(h_textbox2);
    % -- color box
    delete(h_rvd1_region);
    delete(h_rvd2_region);
    delete(h_rvd3_region);
    delete(h_rvd4_region);
    delete(h_rvd5_region);
    delete(h_rvd6_region);
    
end
clearvars i k

close(v);
el_time_video = toc/60;
