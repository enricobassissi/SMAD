%% --------------------------------------------------------------------- %%
%% ----------------------- Plot LegWise Orbit -------------------------- %%
%% ------------------------ Architecture 1RL --------------------------- %%
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
load(imp_path+'Asteroids\Workspaces\ws_1RL_0.7_27mN.mat')

sol.departure_mjd2000 = x(1)*sim.TU/86400;

%% orbit plots
% transfer orbits
r_transf_orbit_1  = [output.r.leg1.*cos(output.theta.leg1), ...
    output.r.leg1.*sin(output.theta.leg1), output.z.leg1];
R_transf_orbit_1 = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_1,sim.n_sol,output.Href.leg1);

% Coasting on Ast 1
start_ast1_sec = (sol.departure_mjd2000+output.tEnd.Leg11)*86400;
end_ast1_sec = (sol.departure_mjd2000+output.tEnd.Leg12)*86400;
time_int_ct1 = linspace(start_ast1_sec,end_ast1_sec,50);
y0_coast_ast1 = [r_encounter.astA1*sim.DU; v_encounter.astA1*sim.DU/sim.TU]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[tC1,yC1] = ode113(@rates, time_int_ct1,y0_coast_ast1,options,'sun');

r_transf_orbit_2  = [output.r.leg2.*cos(output.theta.leg2), ...
    output.r.leg2.*sin(output.theta.leg2), output.z.leg2];
R_transf_orbit_2 = rotate_local2ecplitic(r_encounter.astD1,r_transf_orbit_2,sim.n_sol,output.Href.leg2);

% Coasting on Ast 2
start_ast2_sec = (sol.departure_mjd2000+output.tEnd.Leg21)*86400;
end_ast2_sec = (sol.departure_mjd2000+output.tEnd.Leg22)*86400;
time_int_ct2 = linspace(start_ast2_sec,end_ast2_sec,50);
y0_coast_ast2 = [r_encounter.astA2*sim.DU; v_encounter.astA2*sim.DU/sim.TU]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[tC2,yC2] = ode113(@rates, time_int_ct2,y0_coast_ast2,options,'sun');

r_transf_orbit_3  = [output.r.leg3.*cos(output.theta.leg3), ...
    output.r.leg3.*sin(output.theta.leg3), output.z.leg3];
R_transf_orbit_3 = rotate_local2ecplitic(r_encounter.astD2,r_transf_orbit_3,sim.n_sol,output.Href.leg3);

% Coasting on Ast 3
start_ast3_sec = (sol.departure_mjd2000+output.tEnd.Leg31)*86400;
end_ast3_sec = (sol.departure_mjd2000+output.tEnd.Leg32)*86400;
time_int_ct3 = linspace(start_ast3_sec,end_ast3_sec,50);
y0_coast_ast3 = [r_encounter.astA3*sim.DU; v_encounter.astA3*sim.DU/sim.TU]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[tC3,yC3] = ode113(@rates, time_int_ct3, y0_coast_ast3,options,'sun');

r_transf_orbit_4  = [output.r.leg4.*cos(output.theta.leg4), ...
    output.r.leg4.*sin(output.theta.leg4), output.z.leg4];
R_transf_orbit_4 = rotate_local2ecplitic(r_encounter.astD3,r_transf_orbit_4,sim.n_sol,output.Href.leg4);

%%
figure('Name','Leg1')
plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
    'Color',colors(1,:),'DisplayName','Traj Leg1')
hold on
% hc1 = plot3( yC1(:,1)/sim.DU, yC1(:,2)/sim.DU, yC1(:,3)/sim.DU,'*','Color',colors(2,:),...
%     'Markersize',3,'DisplayName','CT');
% hc1.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot3(r_encounter.astA1(1),r_encounter.astA1(2),r_encounter.astA1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.astD1(1),r_encounter.astD1(2),r_encounter.astD1(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_1)

plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth

fraction_of_the_orbit = 1;
hello_orbit1 = sol.departure_mjd2000 + output.tEnd.Leg12;
plot_asteorid_orbit(hello_orbit1,fraction_of_the_orbit,sol.asteroid_1,colors,3);

% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
legend('Location','northoutside')
legend('NumColumns',3)
view(2)
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]');

%%
figure('Name','Leg2')
hold on
hpt2 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
    'Color',colors(1,:),'DisplayName','Traj Leg2');
% hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';
% hc2 = plot3( yC2(:,1)/sim.DU, yC2(:,2)/sim.DU, yC2(:,3)/sim.DU,'*','Color',colors(2,:),...
%     'Markersize',3,'DisplayName','CT');
% hc2.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot3(r_encounter.astD1(1),r_encounter.astD1(2),r_encounter.astD1(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.astA2(1),r_encounter.astA2(2),r_encounter.astA2(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2)
plot3(r_encounter.astD2(1),r_encounter.astD2(2),r_encounter.astD2(3),...
    '*','Color',colors(4,:),'DisplayName',sol.asteroid_2)

plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth

fraction_of_the_orbit = 1;
hello_orbit2 = sol.departure_mjd2000 + output.tEnd.Leg22;
plot_asteorid_orbit(hello_orbit2,fraction_of_the_orbit,sol.asteroid_2,colors,4);

% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
legend('Location','northoutside')
legend('NumColumns',3)
view(2)
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]');

%%
figure('Name','Leg3')
hold on
hpt3 = plot3(R_transf_orbit_3(:,1),R_transf_orbit_3(:,2),R_transf_orbit_3(:,3),...
    'Color',colors(1,:),'DisplayName','Traj Leg3');
% hpt3.Annotation.LegendInformation.IconDisplayStyle = 'off';
% hc3 = plot3( yC3(:,1)/sim.DU, yC3(:,2)/sim.DU, yC3(:,3)/sim.DU,'*','Color',colors(2,:),...
%     'Markersize',3,'DisplayName','CT');
% hc3.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot3(r_encounter.astD2(1),r_encounter.astD2(2),r_encounter.astD2(3),...
    '*','Color',colors(4,:),'DisplayName',sol.asteroid_2)
plot3(r_encounter.astA3(1),r_encounter.astA3(2),r_encounter.astA3(3),...
    '^','Color',colors(5,:),'DisplayName',sol.asteroid_3)
plot3(r_encounter.astD3(1),r_encounter.astD3(2),r_encounter.astD3(3),...
    '*','Color',colors(5,:),'DisplayName',sol.asteroid_3)

plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth

fraction_of_the_orbit = 1;
hello_orbit3 = sol.departure_mjd2000 + output.tEnd.Leg32;
plot_asteorid_orbit(hello_orbit3,fraction_of_the_orbit,sol.asteroid_3,colors,5);

% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
legend('Location','northoutside')
legend('NumColumns',3)
view(2)
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]');

%%
figure('Name','Leg4')
hold on
hpt4 = plot3(R_transf_orbit_4(:,1),R_transf_orbit_4(:,2),R_transf_orbit_4(:,3),...
    'Color',colors(1,:),'DisplayName','Traj Leg4');
% hpt4.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot3(r_encounter.astD3(1),r_encounter.astD3(2),r_encounter.astD3(3),...
    '*','Color',colors(5,:),'DisplayName',sol.asteroid_3)
plot3(r_encounter.astA4(1),r_encounter.astA4(2),r_encounter.astA4(3),...
    '^','Color',colors(6,:),'DisplayName',sol.asteroid_4)
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
% PLANETS
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
% Asteroids
fraction_of_the_orbit = 1;
hello_orbit4 = sol.departure_mjd2000 + output.tEnd.Leg41;
plot_asteorid_orbit(hello_orbit4,fraction_of_the_orbit,sol.asteroid_4,colors,6);
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
legend('Location','northoutside')
legend('NumColumns',3)
view(2)