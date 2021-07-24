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

%% --------------------------------------------------------------------- %%
%% ------------------------------ SC1 ---------------------------------- %%
%% --------------------------------------------------------------------- %%
%% LOAD INIT GUESS DATA
load('160kg_dry_64mN.mat')
load('Mission_160kg_dry.mat')
load('MA5_160kg_dry.mat')
load('data_before_video_sc1.mat')

%% aspect angles
% -- position
STR_SA_SC1.SC_Traj = [xpl,ypl,zpl]*sim.DU;
% -- velocity
[v_cart_local_SC1_leg1] = from_dt_sol_to_cart(SC1.leg1.HS.X);
[v_cart_local_SC1_leg2] = from_dt_sol_to_cart(SC1.leg2.HS.X);
SC1_Vel = [v_cart_local_SC1_leg1; SC1.coasting.leg1.v_ast; 
    v_cart_local_SC1_leg2; SC1.coasting.leg2.v_ast;]; % AU
SC1.uniform.V_cart = SC1_Vel;

time_vect = SC1.uniform.time;
time_vect(100) = (time_vect(101)+time_vect(99))/2;
time_vect(201) = (time_vect(202)+time_vect(200))/2;
time_vect(299) = time_vect(298)+(time_vect(301)-time_vect(298))/3;
time_vect(300) = time_vect(298)+(time_vect(301)-time_vect(298))*2/3;

dt_cruise = time_vect(2)-time_vect(1);
dt_rv = time_vect(103)-time_vect(102);
number_of_new_time_intervals = round((time_vect(end)-time_vect(1))/dt_rv);
new_time_vect = linspace(time_vect(1),time_vect(end),number_of_new_time_intervals)';

STR_SA_SC1.time = new_time_vect;

[v_az_SC1_sph_pre,v_elev_SC1_sph_pre,v_r_SC1_sph_pre] = cart2sph(SC1_Vel(:,1),SC1_Vel(:,2),SC1_Vel(:,3));
v_az_SC1_sph_post = smooth(interp1(time_vect,...
    smooth(unwrap(v_az_SC1_sph_pre),6),new_time_vect),20);
v_elev_SC1_sph_post = interp1(time_vect,smooth(v_elev_SC1_sph_pre,5),new_time_vect);
v_r_SC1_sph_post = interp1(time_vect,smooth(v_r_SC1_sph_pre,5),new_time_vect);
[vxpl,vypl,vzpl] = sph2cart(wrapToPi(v_az_SC1_sph_post),v_elev_SC1_sph_post,v_r_SC1_sph_post);

STR_SA_SC1.SC_Vel = [vxpl,vypl,vzpl];

[SC1.angles.SAA,SC1.angles.EVA,SC1.angles.SCA,SC1.angles.SolarConjunction,...
    SC1.angles.y_SC_EA] = aspect_angles_LT(STR_SA_SC1);
SC1.angles.time_mjd2000 = STR_SA_SC1.time;

% earth position at same moments, for melia
clearvars rEA vEA kepEA
for i=1:length(SC1.angles.time_mjd2000)
    [kepEA,ksun] = uplanet(SC1.angles.time_mjd2000(i),3);
    [rEA(i,:),vEA(i,:)] = sv_from_coe(kepEA,ksun);
end
SC1.angles.EarthPosition = rEA;
% figure()
% plot3(rEA(:,1),rEA(:,2),rEA(:,3))
% figure()
% plot(SC1.angles.time_mjd2000,rEA)

%% Plot Angles SC1
figure()
plot((STR_SA_SC1.time-sol.departure_mjd2000),rad2deg(SC1.angles.SCA),'Color',colors(1,:))
xlabel('time [d]'); ylabel('Sun Conjunction Angle [deg]');
hold on
% plot(STR_SA_SC1.time,SC1.angles.SolarConjunction)
xlined = [];
for i=1:length(SC1.angles.SolarConjunction)
    if SC1.angles.SolarConjunction(i) == 1
        xlined=[xlined,(STR_SA_SC1.time(i)-sol.departure_mjd2000)];
    end
end
xline(xlined,'Color',colors(6,:))
xlabel('time [d]'); ylabel('Sun Conjunction Angle [deg]');

figure()
plot((STR_SA_SC1.time-sol.departure_mjd2000),rad2deg(SC1.angles.SAA),'Color',colors(1,:))
xlabel('time [d]'); ylabel('Sun Aspect Angle [deg]');

figure()
plot((STR_SA_SC1.time-sol.departure_mjd2000),rad2deg(SC1.angles.EVA),'Color',colors(1,:))
xlabel('time [d]'); ylabel('Earth Visibility Angle [deg]');

figure()
plot((STR_SA_SC1.time-sol.departure_mjd2000),vecnorm(SC1.angles.y_SC_EA,2,2)/sim.DU)
xlabel('time [d]'); ylabel('Rel dist SC EA [AU]');

%% --------------------------------------------------------------------- %%
%% ------------------------------ SC2 ---------------------------------- %%
%% --------------------------------------------------------------------- %%
%% LOAD INIT GUESS DATA
clear; close; clc;
load('160kg_dry_64mN.mat')
load('Mission_160kg_dry.mat')
load('MA5_160kg_dry.mat')
load('data_before_video_sc2.mat')

%% aspect angles
% -- position
STR_SA_SC2.SC_Traj = [xpl,ypl,zpl]*sim.DU;
% -- velocity
[v_cart_local_SC2_lega] = from_dt_sol_to_cart(SC2.lega.HS.X);
[v_cart_local_SC2_legb] = from_dt_sol_to_cart(SC2.legb.HS.X);
SC2_Vel = [v_cart_local_SC2_lega; SC2.coasting.lega.v_ast; 
    v_cart_local_SC2_legb; SC2.coasting.legb.v_ast;]; % AU
SC2.uniform.V_cart = SC2_Vel;

time_vect = SC2.uniform.time;
time_vect(99) = time_vect(98)+(time_vect(101)-time_vect(98))/3;
time_vect(100) = time_vect(98)+(time_vect(101)-time_vect(98))*2/3;
time_vect(201) = (time_vect(202)+time_vect(200))/2;
time_vect(300) = (time_vect(301)+time_vect(299))/2;

dt_cruise = time_vect(2)-time_vect(1);
dt_rv = time_vect(103)-time_vect(102);
number_of_new_time_intervals = round((time_vect(end)-time_vect(1))/dt_rv);
new_time_vect = linspace(time_vect(1),time_vect(end),number_of_new_time_intervals)';

STR_SA_SC2.time = new_time_vect;

[v_az_SC2_sph_pre,v_elev_SC2_sph_pre,v_r_SC2_sph_pre] = cart2sph(SC2_Vel(:,1),SC2_Vel(:,2),SC2_Vel(:,3));
v_az_SC2_sph_post = smooth(interp1(time_vect,smooth(unwrap(v_az_SC2_sph_pre),6),new_time_vect),20);
v_elev_SC2_sph_post = interp1(time_vect,smooth(v_elev_SC2_sph_pre,5),new_time_vect);
v_r_SC2_sph_post = interp1(time_vect,smooth(v_r_SC2_sph_pre,5),new_time_vect);
[vxpl,vypl,vzpl] = sph2cart(wrapToPi(v_az_SC2_sph_post),v_elev_SC2_sph_post,v_r_SC2_sph_post);

STR_SA_SC2.SC_Vel = [vxpl,vypl,vzpl];

[SC2.angles.SAA,SC2.angles.EVA,SC2.angles.SCA,SC2.angles.SolarConjunction,...
    SC2.angles.y_SC_EA] = aspect_angles_LT(STR_SA_SC2);
SC2.angles.time_mjd2000 = STR_SA_SC2.time;

clearvars rEA vEA kepEA
for i=1:length(SC2.angles.time_mjd2000)
    [kepEA,ksun] = uplanet(SC2.angles.time_mjd2000(i),3);
    [rEA(i,:),vEA(i,:)] = sv_from_coe(kepEA,ksun);
end
SC2.angles.EarthPosition = rEA;
% figure()
% plot3(rEA(:,1),rEA(:,2),rEA(:,3))
% figure()
% plot(SC2.angles.time_mjd2000,rEA)

%% Plot Angles SC2
figure()
plot((STR_SA_SC2.time-sol.departure_mjd2000),rad2deg(SC2.angles.SCA),'Color',colors(1,:))
xlabel('time [d]'); ylabel('Sun Conjunction Angle [deg]');
hold on
% plot(sol.SCtime,sol.angles.SolarConjunction)
xlined = [];
for i=1:length(SC2.angles.SolarConjunction)
    if SC2.angles.SolarConjunction(i) == 1
        xlined=[xlined,(STR_SA_SC2.time(i)-sol.departure_mjd2000)];
    end
end
% xline(xlined,'Color',colors(6,:))
xlabel('time [d]'); ylabel('Sun Conjunction Angle [deg]');

figure()
plot((STR_SA_SC2.time-sol.departure_mjd2000),rad2deg(SC2.angles.SAA),'Color',colors(1,:))
xlabel('time [d]'); ylabel('Sun Aspect Angle [deg]');

figure()
plot((STR_SA_SC2.time-sol.departure_mjd2000),rad2deg(SC2.angles.EVA),'Color',colors(1,:))
xlabel('time [d]'); ylabel('Earth Visibility Angle [deg]');

figure()
plot((STR_SA_SC2.time-sol.departure_mjd2000),vecnorm(SC2.angles.y_SC_EA,2,2)/sim.DU)
xlabel('time [d]'); ylabel('Rel dist SC EA [AU]');