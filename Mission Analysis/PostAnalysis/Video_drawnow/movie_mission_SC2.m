%% MOVIE PLOT OF THE WHOLE MISSION SC2
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

%% --------------------------------------------------------------------- %%
%% --------------------------------- SC2 ------------------------------- %%
%% --------------------------------------------------------------------- %%
%% the time vector
% MA.SC1.uniform.time = [SC1.mjd2000_dep+SC1.leg1.timead.*sim.TU/86400; SC1.coasting.leg1.time;
%     SC1.coasting.leg1.time(end)+SC1.leg2.timead.*sim.TU/86400; SC1.coasting.leg2.time;]; % mjd 2000
time_vect = MA.SC2.mjd2000_dep+MA.SC2.uniform.time/86400;
err_1st_tof = time_vect(101)-time_vect(100)
new_tofa = time_vect(100)-time_vect(1)
new_cta = SC1.CT1+time_vect(101)-time_vect(100)
err_2nd_tof = time_vect(301)-time_vect(300)
new_tofb = time_vect(300)-time_vect(201)

% time_vect2 = MA.SC2.mjd2000_dep+MA.SC2.uniform.time/86400;
% time_vect2(101)-time_vect2(100)
% time_vect2(301)-time_vect2(300)
% -- length
n = length(time_vect);

%% sc trajectory
% --- new encounters
MJD01_dim = SC2.mjd2000_dep;
[kep_EA,~] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,sim.mu_dim);
r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;
new_r_encounter.EA = r_EA;
new_v_encounter.EA = v_EA;

MJDAa_dim = (SC2.mjd2000_dep+new_tofa);
[kep_ast_Aa] = uNEO3(MJDAa_dim,SC2.asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[rAa, vAa] = sv_from_coe(kep_ast_Aa,sim.mu_dim); % km, km/s
rAa = rAa/sim.DU;
vAa = vAa/sim.DU*sim.TU;
new_r_encounter.astAa = rAa;
new_v_encounter.astAa = vAa;

MJDDa_dim = (SC2.mjd2000_dep+new_tofa+new_cta);
[kep_ast_Da] = uNEO3(MJDDa_dim,SC2.asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[rDa, vDa] = sv_from_coe(kep_ast_Da,sim.mu_dim); % km, km/s
rDa = rDa/sim.DU;
vDa = vDa/sim.DU*sim.TU;
new_r_encounter.astDa = rDa;
new_v_encounter.astDa = vDa;

MJDAb_dim = (SC2.mjd2000_dep+new_tofa+new_cta+new_tofb);
[kep_ast_Ab] = uNEO3(MJDAb_dim,SC2.asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
[rAb, vAb] = sv_from_coe(kep_ast_Ab,sim.mu_dim); % km, km/s
rAb = rAb/sim.DU;
vAb = vAb/sim.DU*sim.TU;
new_r_encounter.astAb = rAb;
new_v_encounter.astAb = vAb;

MJDDb_dim = (SC2.mjd2000_dep+new_tofa+new_cta+new_tofb+new_cta);
[kep_ast_Db] = uNEO3(MJDDb_dim,SC2.asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
[rDb, vDb] = sv_from_coe(kep_ast_Db,sim.mu_dim); % km, km/s
rDb = rDb/sim.DU;
vDb = vDb/sim.DU*sim.TU;
new_r_encounter.astDb = rDb;
new_v_encounter.astDb = vDb;

%% DT - SC 2 - LEG a
data.Tmax = 0.015; % max thrust available for EGI
data.Pmax = 450; % max power related to 15 mN thrust on qitetiq t5
data.Isp = 2750/sim.TU; % for 15 mN the Isp of Qinetiq T5 is 2750s
data.ThrustModulationFlag = 2;
data.angle_inplane_panels_max = 10; % degrees
data.n_int = 50;
data.muS = sim.mu_dim;   
data.MU = sol.M_start_SC2_lega; % mass adimensionalisation on wet mass
SC2.asteroid_a = sol.asteroid_a;

v_launcher = sol.v_inf_magn/sim.DU*sim.TU*[cos(sol.el)*cos(sol.az); cos(sol.el)*sin(sol.az); sin(sol.el)];
v_dep = new_v_encounter.EA + v_launcher;  %if parabolic escape (v_extra = 0)
SC2.lega.v_in = v_dep;
SC2.lega.v_end = new_v_encounter.astAa;
SC2.lega.r_in = new_r_encounter.EA;
SC2.lega.r_end = new_r_encounter.astAa;
SC2.lega.N_rev = sol.Nrev(3);
SC2.lega.TOF = new_tofa*86400/sim.TU;
SC2.lega.Href = sol.Href(:,3);

DT_input = SC2.lega;
tic
[SC2.lega.HS, SC2.lega.RES, SC2.lega.timead, SC2.lega.Xad, SC2.lega.Xpropad, ...
    SC2.lega.Xpropad_DT] = DT_executable(DT_input, sim, data);
el_time.DT_SC2_lega = toc/60; % minutes

SC2.lega.PROPULSION = nikita(SC2.lega, data, sim, colors);
SC2.lega.EPS = marco_elia(SC2.lega,new_r_encounter.EA,data,sim,SC2.lega.Href,colors);

%% -- if run with less than 100 points
data.n_int = 100;
[SC2.lega] = interpoliamo_porc(SC2.lega,data.n_int);
SC2.lega.v_in = v_dep;
SC2.lega.v_end = new_v_encounter.astAa;
SC2.lega.r_in = new_r_encounter.EA;
SC2.lega.r_end = new_r_encounter.astAa;
SC2.lega.N_rev = sol.Nrev(3);
SC2.lega.TOF = new_tofa*86400/sim.TU;
SC2.lega.Href = sol.Href(:,3);
SC2.mjd2000_dep = sol.departure_mjd2000;

SC2.lega.PROPULSION = nikita(SC2.lega, data, sim, colors);
SC2.lega.EPS = marco_elia(SC2.lega,new_r_encounter.EA,data,sim,SC2.lega.Href,colors);
% -- new mass after DT calculations
SC2.lega.mass_start = SC2.lega.HS.X(1,7);
SC2.lega.mass_end = SC2.lega.HS.X(end,7);
SC2.lega.mass_depleted = SC2.lega.mass_start - SC2.lega.mass_end;
SC2.lega.mass_fraction = SC2.lega.mass_depleted/SC2.lega.mass_start;

% -- coasting between the two legs
[SC2.coasting.lega] = relative_position_coasting_stuff(SC2.asteroid_a,SC2.mjd2000_dep,...
    new_tofa,new_cta,data.n_int,EphData,sim,colors);

%% DT - SC 2 - LEG b
data.Tmax = 0.02; % max thrust available for EGI
data.Pmax = 550; % max power related to 15 mN thrust on qitetiq t5
data.Isp = 2900/sim.TU; % for 15 mN the Isp of Qinetiq T5 is 2750s
data.ThrustModulationFlag = 1;
data.angle_inplane_panels_max = 20; % degrees
data.n_int = 20;
% data.MU = sol.M_start_SC2_legb; % mass adimensionalisation on wet mass
data.MU = SC2.lega.mass_end - sim.M_pods; % mass adimensionalisation on wet mass
SC2.asteroid_b = sol.asteroid_b;

SC2.legb.v_in = new_v_encounter.astDa;
SC2.legb.v_end = new_v_encounter.astAb;
SC2.legb.r_in = new_r_encounter.astDa;
SC2.legb.r_end = new_r_encounter.astAb;
SC2.legb.N_rev = sol.Nrev(4);
SC2.legb.TOF = new_tofb*86400/sim.TU;
SC2.legb.Href = sol.Href(:,4);

DT_input = SC2.legb;
tic
[SC2.legb.HS, SC2.legb.RES, SC2.legb.timead, SC2.legb.Xad, SC2.legb.Xpropad, ...
    SC2.legb.Xpropad_DT] = DT_executable(DT_input, sim, data);
el_time.DT_SC2_legb = toc/60; % minutes

SC2.legb.PROPULSION = nikita(SC2.legb, data, sim, colors);
SC2.legb.EPS = marco_elia(SC2.legb,new_r_encounter.astDa,data,sim,SC2.legb.Href,colors);

%% -- if run with less than 100 points
data.n_int = 100;
[SC2.legb] = interpoliamo_porc(SC2.legb,data.n_int);
SC2.legb.v_in = new_v_encounter.astDa;
SC2.legb.v_end = new_v_encounter.astAb;
SC2.legb.r_in = new_r_encounter.astDa;
SC2.legb.r_end = new_r_encounter.astAb;
SC2.legb.N_rev = sol.Nrev(4);
SC2.legb.TOF = new_tofb*86400/sim.TU;
SC2.legb.Href = sol.Href(:,4);
SC2.mjd2000_dep = sol.departure_mjd2000;

SC2.legb.PROPULSION = nikita(SC2.legb, data, sim, colors);
SC2.legb.EPS = marco_elia(SC2.legb,new_r_encounter.astDa,data,sim,SC2.legb.Href,colors);

% -- new mass after DT calculations
SC2.legb.mass_start = SC2.legb.HS.X(1,7);
SC2.legb.mass_end = SC2.legb.HS.X(end,7);
SC2.legb.mass_depleted = SC2.legb.mass_start - SC2.legb.mass_end;
SC2.legb.mass_fraction = SC2.legb.mass_depleted/SC2.legb.mass_start;

% -- coasting on the last asteroid, same coasting time
SC2.mjd2000_dep_asta = SC2.mjd2000_dep+new_tofa+new_cta;
[SC2.coasting.legb] = relative_position_coasting_stuff(SC2.asteroid_b,SC2.mjd2000_dep_asta,...
    new_tofb,new_cta,data.n_int,EphData,sim,colors);

%% -- building of a uniform position, time, thrust vector for SC1
SC2.uniform.R = [SC2.lega.EPS.R_cartesian; SC2.coasting.lega.r_ast; 
    SC2.legb.EPS.R_cartesian; SC2.coasting.legb.r_ast;]./sim.DU; % AU
figure()
% plot3(SC1.uniform.R(:,1),SC1.uniform.R(:,2),SC1.uniform.R(:,3));
plot3(SC2.uniform.R(1:sim.n_sol,1),SC2.uniform.R(1:sim.n_sol,2),SC2.uniform.R(1:sim.n_sol,3),...
    'Color',colors(1,:),'DisplayName','Cruise - Leg a');
axis equal; grid on; hold on;
plot3(SC2.uniform.R(sim.n_sol+1:2*sim.n_sol,1),SC2.uniform.R(sim.n_sol+1:2*sim.n_sol,2),...
    SC2.uniform.R(sim.n_sol+1:2*sim.n_sol,3),...
    'Color',colors(2,:),'DisplayName','Rendezvous - 2011BP40');
plot3(SC2.uniform.R(2*sim.n_sol+1:3*sim.n_sol,1),SC2.uniform.R(2*sim.n_sol+1:3*sim.n_sol,2),...
    SC2.uniform.R(2*sim.n_sol+1:3*sim.n_sol,3),...
    'Color',colors(3,:),'DisplayName','Cruise - Leg b');
plot3(SC2.uniform.R(3*sim.n_sol+1:4*sim.n_sol,1),SC2.uniform.R(3*sim.n_sol+1:4*sim.n_sol,2),...
    SC2.uniform.R(3*sim.n_sol+1:4*sim.n_sol,3),...
    'Color',colors(4,:),'DisplayName','Rendezvous - 2021JE1');
plot3(new_r_encounter.EA(1), new_r_encounter.EA(2), new_r_encounter.EA(3),'*',...
    'Color',colors(8,:),'DisplayName','Earth Dep');
plot3(new_r_encounter.astAa(1), new_r_encounter.astAa(2), new_r_encounter.astAa(3),'^',...
    'DisplayName','2011BP40 Arr');
plot3(new_r_encounter.astDa(1), new_r_encounter.astDa(2), new_r_encounter.astDa(3),'*',...
    'DisplayName','2011BP40 Dep');
plot3(new_r_encounter.astAb(1), new_r_encounter.astAb(2), new_r_encounter.astAb(3),'^',...
    'DisplayName','2021JE1 Arr');
plot3(new_r_encounter.astDb(1), new_r_encounter.astDb(2), new_r_encounter.astDb(3),'*',...
    'DisplayName','End Mission');
hsun = plot3(0,0,0,'*','Color',colors(4,:));
hsun.Annotation.LegendInformation.IconDisplayStyle = 'off';

[~, h_earth_whole] = plot_object_orbit(new_time_vect(1),'earth',365,sim,data,100,colors,8);
[~, h_mars_whole] = plot_object_orbit(new_time_vect(1),'mars',687,sim,data,100,colors,6);
h_earth_whole.Annotation.LegendInformation.IconDisplayStyle = 'off';
h_mars_whole.Annotation.LegendInformation.IconDisplayStyle = 'off';
h_earth_whole.Color(4) = 0.2;
h_mars_whole.Color(4) = 0.2;

view(2)
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
legend('show','Location','northeastoutside');

SC2.uniform.time = [SC2.mjd2000_dep+SC2.lega.timead.*sim.TU/86400; SC2.coasting.lega.time;
    SC2.coasting.lega.time(end)+SC2.legb.timead.*sim.TU/86400; SC2.coasting.legb.time;]; % mjd 2000

SC2.uniform.Thrust = [SC2.lega.EPS.T_transf_orbit; zeros(sim.n_sol,3);
    SC2.legb.EPS.T_transf_orbit; zeros(sim.n_sol,3)]; % N
figure()
plot(SC2.uniform.time, SC2.uniform.Thrust)
xlabel('t [mjd2000]'); ylabel('T [N]'); legend('Tx','Ty','Tz');
figure()
plot(SC2.uniform.time, vecnorm(SC2.uniform.Thrust,2,2))
xlabel('t [mjd2000]'); ylabel('T magnitude [N]');

%% time interpolation
time_vect = SC2.uniform.time;
time_vect(99) = time_vect(98)+(time_vect(101)-time_vect(98))/3;
time_vect(100) = time_vect(98)+(time_vect(101)-time_vect(98))*2/3;
time_vect(201) = (time_vect(202)+time_vect(200))/2;
time_vect(300) = (time_vect(301)+time_vect(299))/2;

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
    [kep_asta] = uNEO3(new_time_vect(k),MA.SC2.asteroid_a,data);
    [r_asta(k,1:3),v_asta] = sv_from_coe(kep_asta,ksun); % km, km/s
    % --- ast 2
    [kep_astb] = uNEO3(new_time_vect(k),MA.SC2.asteroid_b,data);
    [r_astb(k,1:3),v_astb] = sv_from_coe(kep_astb,ksun); % km, km/s
end
clearvars k

% -- adimensionalise
% rME=rME/sim.DU; rVE=rVE/sim.DU; 
rEA=rEA/sim.DU; rMA=rMA/sim.DU;
r_asta = r_asta/sim.DU; r_astb = r_astb/sim.DU;

%% sc trajectory interpolation
[az_SC2_sph_pre,elev_SC2_sph_pre,r_SC2_sph_pre] = cart2sph(SC2.uniform.R(:,1),SC2.uniform.R(:,2),SC2.uniform.R(:,3));
az_SC2_sph_post = smooth(interp1(time_vect,smooth(unwrap(az_SC2_sph_pre),6),new_time_vect),20);
elev_SC2_sph_post = interp1(time_vect,smooth(elev_SC2_sph_pre,5),new_time_vect);
r_SC2_sph_post = interp1(time_vect,smooth(r_SC2_sph_pre,5),new_time_vect);
[xpl,ypl,zpl] = sph2cart(wrapToPi(az_SC2_sph_post),elev_SC2_sph_post,r_SC2_sph_post);
% yplot = 0;
% yplot = SC1.uniform.R;
yplot = [xpl,ypl,zpl];

figure()
hold on
plot(time_vect,unwrap(az_SC2_sph_pre))
plot(new_time_vect,az_SC2_sph_post)
% plot(time_vect,smooth(unwrap(az_SC1_sph_pre),10))
% plot(time_vect,smooth(unwrap(az_SC1_sph_pre),6))
title('az')

figure()
plot(time_vect,elev_SC2_sph_pre)
hold on
plot(new_time_vect,elev_SC2_sph_post)
title('elev')

figure()
plot(time_vect,r_SC2_sph_pre)
hold on
plot(new_time_vect,r_SC2_sph_post)
title('r')

figure()
plot(time_vect,SC2.uniform.R(:,1))
hold on
plot(new_time_vect,xpl)
title('x')

figure()
plot(time_vect,SC2.uniform.R(:,2))
hold on
plot(new_time_vect,ypl)
title('y')

figure()
plot(time_vect,SC2.uniform.R(:,3))
hold on
plot(new_time_vect,zpl)
title('z')

%% Thrust stuff
T_SC2 = [SC2.lega.HS.T;zeros(sim.n_sol,1);SC2.legb.HS.T;zeros(sim.n_sol,1)];
alpha_SC2 = [SC2.lega.HS.alpha;zeros(sim.n_sol,1);SC2.legb.HS.alpha;zeros(sim.n_sol,1)];
beta_SC2 = [SC2.lega.HS.beta;zeros(sim.n_sol,1);SC2.legb.HS.beta;zeros(sim.n_sol,1)];

figure()
plot(time_vect,T_SC2)
xlabel('t [mjd2000]'); ylabel('T [N]'); 
figure()
plot(time_vect,alpha_SC2)
xlabel('t [mjd2000]'); ylabel('T [N]'); 
figure()
plot(time_vect,beta_SC2)
xlabel('t [mjd2000]'); ylabel('T [N]'); 

% --- reinterpolate
T_SC2_post = interp1(time_vect,T_SC2,new_time_vect);
alpha_SC2_post = interp1(time_vect,alpha_SC2,new_time_vect);
beta_SC2_post = interp1(time_vect,beta_SC2,new_time_vect);

figure()
plot(time_vect,T_SC2)
hold on
plot(new_time_vect,T_SC2_post)
xlabel('t [mjd2000]'); ylabel('T [N]');

figure()
plot(time_vect,beta_SC2)
hold on
plot(new_time_vect,beta_SC2_post)
xlabel('t [mjd2000]'); ylabel('beta[rad]'); 
figure()
plot(time_vect,alpha_SC2)
hold on
plot(new_time_vect,alpha_SC2_post)
xlabel('t [mjd2000]'); ylabel('alpha[rad]'); 

%% video stuff SC2
load('data_before_video_sc2.mat')

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
figure('Name','Mission Video SC2'); %,'WindowState','maximized'
set(gca,'nextplot','replacechildren');
set(gcf, 'Position', get(0, 'Screensize'));
v = VideoWriter('mission_SC2_2','MPEG-4');
% v = VideoWriter('mission4', 'Uncompressed AVI');
v.Quality = 100;
open(v);  

% set(gca, 'CameraPosition', [0 lim_max lim_max]);
for i=1:13:n % i=1:4:n % i=1:364:n
    % ---- start orbit subplot
    subplot(4,3,[1,2,4,5,7,8,10,11])
    % --- whole orbits
    axis equal; grid on; hold on;
    xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
    [~, h_ast1_whole] = plot_object_orbit(new_time_vect(1),'asteroid',1.5*365,sim,data,100,colors,3,MA.SC2.asteroid_a);
    [~, h_ast2_whole] = plot_object_orbit(new_time_vect(1),'asteroid',1.5*365,sim,data,100,colors,4,MA.SC2.asteroid_b);
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
    [~, h_ast1_highlight] = plot_object_orbit(new_time_vect(i),'asteroid',0.25*365,sim,data,30,colors,3,MA.SC2.asteroid_a);
    [~, h_ast2_highlight] = plot_object_orbit(new_time_vect(i),'asteroid',0.25*365,sim,data,30,colors,4,MA.SC2.asteroid_b);
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
    
    [~, h_ast1_highlight2] = plot_object_orbit(new_time_vect(i),'asteroid',0.08*365,sim,data,10,colors,3,MA.SC2.asteroid_a);
    [~, h_ast2_highlight2] = plot_object_orbit(new_time_vect(i),'asteroid',0.08*365,sim,data,10,colors,4,MA.SC2.asteroid_b);
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
    hast1 = plot3( r_asta(i,1), r_asta(i,2), r_asta(i,3),'*','Color',colors(3,:),...
        'MarkerSize',5,'DisplayName',MA.SC2.asteroid_a);
    hast2 = plot3( r_astb(i,1), r_astb(i,2), r_astb(i,3),'*','Color',colors(4,:),...
        'MarkerSize',5,'DisplayName',MA.SC2.asteroid_b);

    legend('show','Location','northwest')
    
    date_view=mjd20002date(new_time_vect(i));
    title (sprintf('Mission Date %d %d %d', date_view(1:3)));
%     set(get(gca,'title'),'Position',[15 13 1.00011])

    view([lim_max_helio lim_max_helio lim_max_helio])  
    xlim([-lim_max_helio lim_max_helio])
    ylim([-lim_max_helio lim_max_helio])
    zlim([-lim_max_helio lim_max_helio])

    % -- textbox inside the orbit plots
    if new_time_vect(i) < MJDAa_dim
        cs1 = '\begin{tabular}{l}  \textbf{Phase} \\ Cruise - Leg 1 \\ Rendezvous - 2011BP40 \\ Cruise - Leg 2 \\ Rendezvous - 2021JE1 \end{tabular}';
        cs2 = {'\textbf{Mission Time [d]}'; ['T+', num2str(round(new_time_vect(i)-SC2.mjd2000_dep))] ;
            ['T-', num2str(round(MJDAa_dim-new_time_vect(i)))];
            '-'; '-'};
        h_textbox1 = annotation('textbox', [0.28, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        h_textbox2 = annotation('textbox', [0.40, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        set(h_textbox1, 'String', cs1);
        set(h_textbox2, 'String', cs2);
    elseif new_time_vect(i) >= MJDAa_dim && new_time_vect(i) < MJDDa_dim
        cs1 = '\begin{tabular}{l}  \textbf{Phase} \\ Cruise - Leg 1 \\ Rendezvous - 2011BP40 \\ Cruise - Leg 2 \\ Rendezvous - 2021JE1 \end{tabular}';
        cs2 = {'\textbf{Mission Time [d]}'; 'Completed' ; 
            ['T+', num2str(round(new_time_vect(i)-MJDAa_dim))];
            ['T-', num2str(round(MJDDa_dim-new_time_vect(i)))]; '-'};
        h_textbox1 = annotation('textbox', [0.28, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        h_textbox2 = annotation('textbox', [0.40, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        set(h_textbox1, 'String', cs1);
        set(h_textbox2, 'String', cs2);
    elseif new_time_vect(i) >= MJDDa_dim && new_time_vect(i) < MJDAb_dim
        cs1 = '\begin{tabular}{l}  \textbf{Phase} \\ Cruise - Leg 1 \\ Rendezvous - 2011BP40 \\ Cruise - Leg 2 \\ Rendezvous - 2021JE1 \end{tabular}';
        cs2 = {'\textbf{Mission Time [d]}'; 'Completed' ; 'Completed';
            ['T+', num2str(round(new_time_vect(i)-MJDDa_dim))]; 
            ['T-', num2str(round(MJDAb_dim-new_time_vect(i)))]};
        h_textbox1 = annotation('textbox', [0.28, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        h_textbox2 = annotation('textbox', [0.40, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        set(h_textbox1, 'String', cs1);
        set(h_textbox2, 'String', cs2);
	elseif new_time_vect(i) >= MJDAb_dim
        cs1 = '\begin{tabular}{l}  \textbf{Phase} \\ Cruise - Leg 1 \\ Rendezvous - 2011BP40 \\ Cruise - Leg 2 \\ Rendezvous - 2021JE1 \end{tabular}';
        cs2 = {'\textbf{Mission Time [d]}'; 'Completed' ; 'Completed'; 'Completed';
             ['T+', num2str(round(new_time_vect(i)-MJDAb_dim))]};
        h_textbox1 = annotation('textbox', [0.28, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        h_textbox2 = annotation('textbox', [0.40, 0.28, 0.1, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
        set(h_textbox1, 'String', cs1);
        set(h_textbox2, 'String', cs2);
    end
    
    drawnow

    % ------ thrust profile subplot
    subplot(4,3,[3,6])
    h_plot_subplot3 = plot(new_time_vect-SC2.mjd2000_dep,T_SC2_post,'Color',colors(1,:));
    hold on; grid on;
    h_plot_subplot3_point = plot(new_time_vect(i)-SC1.mjd2000_dep,T_SC2_post(i),'o','Color',colors(2,:));
%     xlabel('T+dep [d]'); 
    ylabel('T [N]'); 
    xlim([min(new_time_vect-SC2.mjd2000_dep) max(new_time_vect-SC2.mjd2000_dep)])
    xticks([0 500 1000 1500 2000])
    ylim([min(T_SC2_post) max(T_SC2_post)])
    
    h_rvd1_region = fill([MJDAa_dim-SC2.mjd2000_dep,MJDDa_dim-SC2.mjd2000_dep,...
        MJDDa_dim-SC2.mjd2000_dep,MJDAa_dim-SC2.mjd2000_dep],...
        [min(T_SC2_post),min(T_SC2_post),...
        max(T_SC2_post),max(T_SC2_post)],'red');
    h_rvd1_region.FaceAlpha=0.3;
    h_rvd1_region.FaceColor=colors(10,:);
    h_rvd1_region.EdgeColor='none';
    
    h_rvd2_region = fill([MJDAb_dim-SC2.mjd2000_dep,MJDDb_dim-SC2.mjd2000_dep,...
        MJDDb_dim-SC2.mjd2000_dep,MJDAb_dim-SC2.mjd2000_dep],...
        [min(T_SC2_post),min(T_SC2_post),...
        max(T_SC2_post),max(T_SC2_post)],'red');
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
    h_plot_subplot61 = plot(new_time_vect-SC2.mjd2000_dep,wrapTo180(rad2deg(alpha_SC2_post)),'Color',colors(1,:));
    hold on; grid on;
    h_plot_subplot61_point = plot(new_time_vect(i)-SC2.mjd2000_dep,wrapTo180(rad2deg(alpha_SC2_post(i))),'o','Color',colors(2,:));
%     xlabel('T+dep [d]'); 
    ylabel('Azimuth [deg]'); 
    xlim([min(new_time_vect-SC2.mjd2000_dep) max(new_time_vect-SC2.mjd2000_dep)])
    xticks([0 500 1000 1500 2000])
    ylim([min(rad2deg(alpha_SC2_post)) max(rad2deg(alpha_SC2_post))])
    
    h_rvd3_region = fill([MJDAa_dim-SC2.mjd2000_dep,MJDDa_dim-SC2.mjd2000_dep,...
        MJDDa_dim-SC2.mjd2000_dep,MJDAa_dim-SC2.mjd2000_dep],...
        [min(rad2deg(alpha_SC2_post)),min(rad2deg(alpha_SC2_post)),...
        max(rad2deg(alpha_SC2_post)),max(rad2deg(alpha_SC2_post))],'red');
    h_rvd3_region.FaceAlpha=0.3;
    h_rvd3_region.FaceColor=colors(10,:);
    h_rvd3_region.EdgeColor='none';
    
    h_rvd4_region = fill([MJDAb_dim-SC2.mjd2000_dep,MJDDb_dim-SC2.mjd2000_dep,...
        MJDDb_dim-SC2.mjd2000_dep,MJDAb_dim-SC2.mjd2000_dep],...
        [min(rad2deg(alpha_SC2_post)),min(rad2deg(alpha_SC2_post)),...
        max(rad2deg(alpha_SC2_post)),max(rad2deg(alpha_SC2_post))],'red');
    h_rvd4_region.FaceAlpha=0.3;
    h_rvd4_region.FaceColor=colors(10,:);
    h_rvd4_region.EdgeColor='none';
    
    h_plot_subplot61.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_plot_subplot61_point.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    subplot(4,3,12)
    h_plot_subplot62 = plot(new_time_vect-SC2.mjd2000_dep,wrapTo180(rad2deg(beta_SC2_post)),'Color',colors(1,:));
    hold on; grid on;
    h_plot_subplot62_point = plot(new_time_vect(i)-SC2.mjd2000_dep,wrapTo180(rad2deg(beta_SC2_post(i))),'o','Color',colors(2,:));
    xlabel('T+dep [d]'); ylabel('Elevation [deg]'); 
    xlim([min(new_time_vect-SC2.mjd2000_dep) max(new_time_vect-SC2.mjd2000_dep)])
    xticks([0 500 1000 1500 2000])
    ylim([min(rad2deg(beta_SC2_post)) max(rad2deg(beta_SC2_post))])
    
    h_rvd5_region = fill([MJDAa_dim-SC2.mjd2000_dep,MJDDa_dim-SC2.mjd2000_dep,...
        MJDDa_dim-SC2.mjd2000_dep,MJDAa_dim-SC2.mjd2000_dep],...
        [min(rad2deg(beta_SC2_post)),min(rad2deg(beta_SC2_post)),...
        max(rad2deg(beta_SC2_post)),max(rad2deg(beta_SC2_post))],'red');
    h_rvd5_region.FaceAlpha=0.3;
    h_rvd5_region.FaceColor=colors(10,:);
    h_rvd5_region.EdgeColor='none';
    
    h_rvd6_region = fill([MJDAb_dim-SC2.mjd2000_dep,MJDDb_dim-SC2.mjd2000_dep,...
        MJDDb_dim-SC2.mjd2000_dep,MJDAb_dim-SC2.mjd2000_dep],...
        [min(rad2deg(beta_SC2_post)),min(rad2deg(beta_SC2_post)),...
        max(rad2deg(beta_SC2_post)),max(rad2deg(beta_SC2_post))],'red');
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
