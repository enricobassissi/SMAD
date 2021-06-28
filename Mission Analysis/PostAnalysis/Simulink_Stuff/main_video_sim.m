%% --------------------------------------------------------------------- %%
%% ----------------------- video of the orbits ------------------------- %%
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
str_path=split(pwd, 'PostAnalysis\Simulink_Stuff');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
% str_path=split(pwd, 'Simulink_Stuff');
% imp_path=string(str_path(1));
addpath(genpath(pwd));

%% Call to NASA JPL Horizons to get Asteroid's Ephemerides
% Import module of Python
try 
    module = py.importlib.import_module('neo_api_function');
catch
    copyfile(py_path+'\neo_api_function.py', pwd, 'f'); 
    module = py.importlib.import_module('neo_api_function');
end
clc

%% Loads
load('MA.mat')
load('data_elements_matrix_44_63_2SC')
asteroid_1 = MA.SC1.asteroid_1;
asteroid_2 = MA.SC1.asteroid_2;
asteroid_a = MA.SC2.asteroid_a;
asteroid_b = MA.SC2.asteroid_b;

%% Parameters
ID_planet = 3; % Earth
mjd2000_dep = MA.SC1.mjd2000_dep;
adim.mu_dim    = 132712440018              ; % actractor parameter [km^3 s^-2]
adim.DU = astroConstants(2); % [km]
adim.TU        = (adim.DU^3/adim.mu_dim )^0.5; % time unit [s]
adim.mu        = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
adim.g0 = 9.81*(adim.TU^2/(1000*adim.DU)); % non-dimensional g0

% --- produce mega table with planets and asteroids ephemerides
time_mjd2000_mega_table_eph = mjd2000_dep + linspace(0,6*365,6*365);
for i = 1:length(time_mjd2000_mega_table_eph)
    [kep_MA, muSun] = uplanet(time_mjd2000_mega_table_eph(i),4);
    [planet_eph.MA(i,1:3), planet_eph.MA(i,4:6)] = sv_from_coe(kep_MA, muSun); % km, km/s
    planet_eph.MA(i,7) = time_mjd2000_mega_table_eph(i);
    
    [kep_EA, muSun] = uplanet(time_mjd2000_mega_table_eph(i),3);
    [planet_eph.EA(i,1:3), planet_eph.EA(i,4:6)] = sv_from_coe(kep_EA, muSun); % km, km/s
    planet_eph.EA(i,7) = time_mjd2000_mega_table_eph(i);
    
    [kep_VE, muSun] = uplanet(time_mjd2000_mega_table_eph(i),2);
    [planet_eph.VE(i,1:3), planet_eph.VE(i,4:6)] = sv_from_coe(kep_VE, muSun); % km, km/s
    planet_eph.VE(i,7) = time_mjd2000_mega_table_eph(i);
    
    [kep_ME, muSun] = uplanet(time_mjd2000_mega_table_eph(i),1);
    [planet_eph.ME(i,1:3), planet_eph.ME(i,4:6)] = sv_from_coe(kep_ME, muSun); % km, km/s
    planet_eph.ME(i,7) = time_mjd2000_mega_table_eph(i);
    
    [kep_ast_A1] = uNEO3(time_mjd2000_mega_table_eph(i),asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
    [planet_eph.A1(i,1:3), planet_eph.A1(i,4:6)] = sv_from_coe(kep_ast_A1,muSun); % km, km/s
    planet_eph.A1(i,7) = time_mjd2000_mega_table_eph(i);

    [kep_ast_A2] = uNEO3(time_mjd2000_mega_table_eph(i),asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
    [planet_eph.A2(i,1:3), planet_eph.A2(i,4:6)] = sv_from_coe(kep_ast_A2,muSun); % km, km/s
    planet_eph.A2(i,7) = time_mjd2000_mega_table_eph(i);
    planet_eph.A2(i,7) = time_mjd2000_mega_table_eph(i);
    
    [kep_ast_Aa] = uNEO3(time_mjd2000_mega_table_eph(i),asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
    [planet_eph.Aa(i,1:3), planet_eph.Aa(i,4:6)] = sv_from_coe(kep_ast_Aa,muSun); % km, km/s
    planet_eph.Aa(i,7) = time_mjd2000_mega_table_eph(i);

    [kep_ast_Ab] = uNEO3(time_mjd2000_mega_table_eph(i),asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
    [planet_eph.Ab(i,1:3), planet_eph.Ab(i,4:6)] = sv_from_coe(kep_ast_Ab,muSun); % km, km/s
    planet_eph.Ab(i,7) = time_mjd2000_mega_table_eph(i);

end

% --- mean distance of planets from the sun
MarsDist = 1.524; % AU
EarthDist = 1; % [AU]
VenusDist = 0.7; % [AU]
MercuryDist = 0.387; % [AU]

% -- whole mission
% % --- increment points of spacecraft trajectory
% resample_factor = 3;
% support_MA_time_1 = MA.SC1.coasting.leg1.time;
% support_MA_time_1(1) = (support_MA_time_1(2) + support_MA_time_1(1))/2;
% support_MA_time_2 = MA.SC1.leg2.timead;
% support_MA_time_2(1) = (support_MA_time_2(2) + support_MA_time_2(1))/2;
% support_MA_time_3 = MA.SC1.coasting.leg2.time;
% support_MA_time_3(1) = (support_MA_time_3(2) + support_MA_time_3(1))/2;
% MA_time_correct = [MA.SC1.leg1.timead;MA.SC1.leg1.timead(end)+support_MA_time_1;
%     MA.SC1.leg1.timead(end)+MA.SC1.coasting.leg1.time(end)+support_MA_time_2;
%     MA.SC1.leg1.timead(end)+MA.SC1.coasting.leg1.time(end)+MA.SC1.leg2.timead(end)+support_MA_time_3];
% % time_mjd2000_traj = mjd2000_dep + linspace(MA_time_correct(1),MA_time_correct(end),length(MA_time_correct)*resample_factor)/86400; % mjd2000, days of the mission
% time_mjd2000_traj = mjd2000_dep + linspace(MA_time_correct(1),MA_time_correct(end),length(MA_time_correct)*resample_factor)/86400; % mjd2000, days of the mission
% % TrajSC_x = interp(MA.SC1.uniform.R(:,1),resample_factor);
% % TrajSC_y = interp(MA.SC1.uniform.R(:,2),resample_factor);
% % TrajSC_z = interp(MA.SC1.uniform.R(:,3),resample_factor);
% TrajSC = interp1(mjd2000_dep + MA_time_correct/86400,MA.SC1.uniform.R,time_mjd2000_traj);
% TrajSC = [TrajSC, time_mjd2000_traj'];
% 
% sim_time = MA_time_correct(end)/86400; % [s]
%
% % -- check plots
% figure('Name','interp orbit check')
% plot3(TrajSC(:,1),TrajSC(:,2),TrajSC(:,3),...
%     'Color',colors(1,:),'DisplayName','Interp')
% hold on; grid on; axis equal;
% plot3(MA.SC1.uniform.R(:,1),MA.SC1.uniform.R(:,2),MA.SC1.uniform.R(:,3),...
%     '--','Color',colors(2,:),'DisplayName','Original')
% xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); 
% legend('show')

% ---- first leg
% --- increment points of spacecraft trajectory
resample_factor = 3;
support_MA_time_1 = MA.SC1.coasting.leg1.time;
support_MA_time_1(1) = (support_MA_time_1(2) + support_MA_time_1(1))/2;
MA_time_correct = [MA.SC1.leg1.timead;MA.SC1.leg1.timead(end)+support_MA_time_1];
% time_mjd2000_traj = mjd2000_dep + linspace(MA_time_correct(1),MA_time_correct(end),length(MA_time_correct)*resample_factor)/86400; % mjd2000, days of the mission
time_mjd2000_traj = mjd2000_dep + linspace(MA_time_correct(1),MA_time_correct(end),length(MA_time_correct)*resample_factor)/86400; % mjd2000, days of the mission
% TrajSC_x = interp(MA.SC1.uniform.R(:,1),resample_factor);
% TrajSC_y = interp(MA.SC1.uniform.R(:,2),resample_factor);
% TrajSC_z = interp(MA.SC1.uniform.R(:,3),resample_factor);
TrajSC = interp1(mjd2000_dep + MA_time_correct/86400,MA.SC1.uniform.R(1:200,:),time_mjd2000_traj);
TrajSC = [TrajSC, time_mjd2000_traj'];

sim_time = MA_time_correct(end)/86400; % [s]

% -- check plots
figure('Name','interp orbit check')
plot3(TrajSC(:,1),TrajSC(:,2),TrajSC(:,3),...
    'Color',colors(1,:),'DisplayName','Interp')
hold on; grid on; axis equal;
plot3(MA.SC1.uniform.R(1:200,1),MA.SC1.uniform.R(1:200,2),MA.SC1.uniform.R(1:200,3),...
    '--','Color',colors(2,:),'DisplayName','Original')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); 
legend('show')

%% Simulation
sim_time = MA.SC1.uniform.time(end)/86400; % [s]
abs_tol = 1e-7;
rel_tol = 1e-12;
solver_name = 'ode113';

load_system('trajectory_plots.slx')
set_param('trajectory_plots','Solver', solver_name,'RelTol', num2str(rel_tol),...
    'AbsTol',num2str(abs_tol),'StopTime',num2str(sim_time));
simulation = sim('trajectory_plots','SimulationMode','accelerator','SrcWorkspace','current');

el_time = toc;
fprintf('Elapsed time: %f min... \n', el_time/60)

% -------------  Extract data from simulink model  -------------------- %
% simulation time
tout = get(simulation,'tout');