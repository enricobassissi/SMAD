%% --------------------------------------------------------------------- %%
%% ----------------------- Plot LegWise Orbit -------------------------- %%
%% ------------------------ Architecture 2RL --------------------------- %%
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
          0    0    0]./255; % (12) BLACK
      
%% add path of functions and python stuff
str_path=split(pwd, 'TrajOptimisation\LowThrust\Asteroids\PlotFunctions\LegWise');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'Asteroids\PlotFunctions\LegWise');
imp_path=string(str_path(1));
addpath(genpath(imp_path));

%% load
load(imp_path+'Asteroids\Ale Workspaces\case22.mat')

sol.departure_mjd2000 = x(1)*sim.TU/86400;

%% orbit plots
% transfer orbits
% SC1
r_transf_orbit_1  = [output.r.leg1.*cos(output.theta.leg1), ...
    output.r.leg1.*sin(output.theta.leg1), output.z.leg1];
R_transf_orbit_1 = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_1,sim.n_sol,output.Href.leg1);

% % Coasting on Ast 1
% start_ast1_sec = (sol.departure_mjd2000+sol.TOF1)*86400;
% end_ast1_sec = (sol.departure_mjd2000+sol.TOF1+sol.CT1)*86400;
% time_int_SC1 = linspace(start_ast1_sec,end_ast1_sec,50);
% y0_coast_ast1 = [r_encounter.astA1*sim.DU; v_encounter.astA1*sim.DU/sim.TU]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
% options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
% [tC1,yC1] = ode113(@rates, time_int_SC1,y0_coast_ast1,options,'sun');

r_transf_orbit_2  = [output.r.leg2.*cos(output.theta.leg2), ...
    output.r.leg2.*sin(output.theta.leg2), output.z.leg2];
R_transf_orbit_2 = rotate_local2ecplitic(r_encounter.astD1,r_transf_orbit_2,sim.n_sol,output.Href.leg2);

% SC2
r_transf_orbit_a  = [output.r.lega.*cos(output.theta.lega), ...
    output.r.lega.*sin(output.theta.lega), output.z.lega];
R_transf_orbit_a = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_a,sim.n_sol,output.Href.lega);

% % Coasting on Ast a
% start_asta_sec = (sol.departure_mjd2000+sol.TOFa)*86400;
% end_asta_sec = (sol.departure_mjd2000+sol.TOFa+sol.CTa)*86400;
% % time_int_SC2 = linspace(start_asta_sec,end_asta_sec,50);
% time_int_SC2 = sol.departure_mjd2000*86400+output.CTa*sim.TU;
% y0_coast_asta = [r_encounter.astAa*sim.DU; v_encounter.astAa*sim.DU/sim.TU]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
% options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
% [tCa,yCa] = ode113(@rates, time_int_SC2, y0_coast_asta,options,'sun');

r_transf_orbit_b  = [output.r.legb.*cos(output.theta.legb), ...
    output.r.legb.*sin(output.theta.legb), output.z.legb];
R_transf_orbit_b = rotate_local2ecplitic(r_encounter.astDa,r_transf_orbit_b,sim.n_sol,output.Href.legb);

%% SC1 - Leg 1
figure('Name','Leg1')
hold on
plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
    'Color',colors(1,:),'DisplayName','Traj SC1 Leg1')
% hc1 = plot3( yC1(:,1)/sim.DU, yC1(:,2)/sim.DU, yC1(:,3)/sim.DU,'*','Color',colors(1,:),...
%     'Markersize',3,'DisplayName','CT');
% hc1.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Earth Dep')
plot3(r_encounter.astA1(1),r_encounter.astA1(2),r_encounter.astA1(3),...
    '^','Color',colors(2,:),'DisplayName',sol.asteroid_1+' Arr')
plot3(r_encounter.astD1(1),r_encounter.astD1(2),r_encounter.astD1(3),...
    '*','Color',colors(2,:),'DisplayName',sol.asteroid_1+' Dep')

axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]');

% PLANETS
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
% plot_planet_orbit(x(1)*sim.TU/(3600*24),4,colors,6); % mars

% Asteroids
fraction_of_the_orbit = 1;
hello_orbit1 = sol.departure_mjd2000+output.t1(end)*sim.TU/(3600*24);
plot_asteorid_orbit(hello_orbit1,fraction_of_the_orbit,sol.asteroid_1,colors,2);

% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
legend('Location','northoutside')
legend('NumColumns',3)
view(2)

%% SC1 - Leg 2
figure('Name','Leg2')
hold on
hpt2 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
    'Color',colors(1,:),'DisplayName','Traj SC1 Leg2');

plot3(r_encounter.astA1(1),r_encounter.astA1(2),r_encounter.astA1(3),...
    '^','Color',colors(2,:),'DisplayName',sol.asteroid_1+' Arr')
plot3(r_encounter.astD1(1),r_encounter.astD1(2),r_encounter.astD1(3),...
    '*','Color',colors(2,:),'DisplayName',sol.asteroid_1+' Dep')
plot3(r_encounter.astA2(1),r_encounter.astA2(2),r_encounter.astA2(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2+' Arr')

axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]');

% PLANETS
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
% plot_planet_orbit(x(1)*sim.TU/(3600*24),4,colors,6); % mars

% Asteroids
fraction_of_the_orbit = 1;
hello_orbit1 = sol.departure_mjd2000+output.t1(end)*sim.TU/(3600*24);
plot_asteorid_orbit(hello_orbit1,fraction_of_the_orbit,sol.asteroid_1,colors,3);
hello_orbit2 = sol.departure_mjd2000+(output.t1(end)+ sol.CT1 + output.t2(end))*sim.TU/(3600*24);
plot_asteorid_orbit(hello_orbit2,fraction_of_the_orbit,sol.asteroid_2,colors,4);

% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
legend('Location','northoutside')
legend('NumColumns',3)
view(2)

%% SC2 - Leg a
figure('Name','Lega')
hold on
plot3(R_transf_orbit_a(:,1),R_transf_orbit_a(:,2),R_transf_orbit_a(:,3),...
    'Color',colors(2,:),'DisplayName','Traj SC2 Leg1');
% hca = plot3( yCa(:,1)/sim.DU, yCa(:,2)/sim.DU, yCa(:,3)/sim.DU,'*','Color',colors(2,:),...
%     'Markersize',3,'DisplayName','CT2');
% hca.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Earth Dep')
plot3(r_encounter.astAa(1),r_encounter.astAa(2),r_encounter.astAa(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_a+' Arr')
plot3(r_encounter.astDa(1),r_encounter.astDa(2),r_encounter.astDa(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_a+' Dep')

axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]');

% PLANETS
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
% plot_planet_orbit(x(1)*sim.TU/(3600*24),4,colors,6); % mars

% Asteroids
fraction_of_the_orbit = 1;
hello_orbita = sol.departure_mjd2000+output.ta(end)*sim.TU/(3600*24);
plot_asteorid_orbit(hello_orbita,fraction_of_the_orbit,sol.asteroid_a,colors,3);

% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
legend('Location','northoutside')
legend('NumColumns',3)
view(2)

%% SC2 - Leg b
figure('Name','Legb')
hold on
hpt4 = plot3(R_transf_orbit_b(:,1),R_transf_orbit_b(:,2),R_transf_orbit_b(:,3),...
    'Color',colors(2,:),'DisplayName','Traj SC2 Leg2');

plot3(r_encounter.astAa(1),r_encounter.astAa(2),r_encounter.astAa(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_a+' Arr')
plot3(r_encounter.astDa(1),r_encounter.astDa(2),r_encounter.astDa(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_a+' Dep')
plot3(r_encounter.astAb(1),r_encounter.astAb(2),r_encounter.astAb(3),...
    '^','Color',colors(9,:),'DisplayName',sol.asteroid_b+' Arr')

axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]');

% PLANETS
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
% plot_planet_orbit(x(1)*sim.TU/(3600*24),4,colors,6); % mars

% Asteroids
fraction_of_the_orbit = 1;
hello_orbita = sol.departure_mjd2000+output.ta(end)*sim.TU/(3600*24);
plot_asteorid_orbit(hello_orbita,fraction_of_the_orbit,sol.asteroid_a,colors,3);
hello_orbitb = sol.departure_mjd2000+(output.ta(end)+ sol.CTa + output.tb(end))*sim.TU/(3600*24);
plot_asteorid_orbit(hello_orbitb,fraction_of_the_orbit,sol.asteroid_b,colors,9);

% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
legend('Location','northoutside')
legend('NumColumns',3)
view(2)
