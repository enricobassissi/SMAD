%% --------------------------------------------------------------------- %%
%% ----------------------- Fixed Adim Points Transfer ------------------ %%
%% ------------------------- ARCH 1+4, LT FLYBY ------------------------ %%
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

sim.case_name = 'ARCH ID 6: LOW THRUST FLYBY ON EVERY ASTEROID';

%% add path of functions and python stuff
str_path=split(pwd, 'TrajOptimisation\LowThrust\JP_LT');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'JP_LT');
imp_path=string(str_path(1));
addpath(genpath(imp_path));

%% simulation parameters
sim.mu_dim    = 132712440018          ; % actractor parameter [km^3 s^-2] -- sun
sim.DU    = 149597870.7           ; % distance unit [km]
sim.TU    = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu    = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol = 1000;                    % number of computational nodes %% 100 initially
sim.x = linspace(0,1,sim.n_sol)';   % 
sim.out_shape = 2;                  % out-of-plane shape (2:CONWAY)
sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = 1;                     % direction of integration (1 FW, -1 BW)
sim.tol_vers = 1e-8;

%% Propulsive system parameters
PS.Is = 3000/sim.TU;  % non-dimensional specific impulse

%% Orbits parameters
RI = [ 1  0 0.1]'; % initial position [DU]  %% ORIGINALE: [1 0 0.1] 
RF = [-1 -1 0]'; % final position [DU]

VI = [ 0  1 0]'; % initial velocity [DU/TU]
VF = [ 0.7  -0.7 0]'; % final velocity [DU/TU]

N_rev = 1; % number of revolution
TOF = 25.4; % TOF [TU] %25.4
M = 1000; % SC mass [kg]
hp = 3; %3 OUT OF PLANE SHAPE PARAMETERS
kp = 3; %3

%% conway

sim.TOF_imposed_flag = 1;

[output] = CW_LowLambert( RI , RF , VI , VF , N_rev , TOF ,M ,hp , kp , PS ,sim );

output.m(1);

figure()
subplot(5,1,1)
plot(output.t*sim.TU/86400,output.u(:,1));
xlabel('Time [days]')
ylabel('In-plane Thrust [N]')

subplot(5,1,2)
plot(output.t*sim.TU/86400,180/pi*output.u(:,2));
xlabel('Time [days]')
ylabel('In-plane Thrust angle [deg]')

subplot(5,1,3)
plot(output.t*sim.TU/86400,output.u(:,3));
xlabel('Time [days]')
ylabel('out-of-plane Thrust [N]')

subplot(5,1,4)
plot(output.t*sim.TU/86400,sqrt(output.u(:,1).^2 + output.u(:,3).^2));
xlabel('Time [days]')
ylabel('Thrust [N]')

subplot(5,1,5)
plot(output.t*sim.TU/86400,output.m);
xlabel('Time [days]')
ylabel('Mass [kg]')

R_sdrp  = [output.r.*cos(output.l) output.r.*sin(output.l) output.z];
Rglobal = rotate_local2ecplitic(RI,R_sdrp,sim.n_sol,output.h_ref_v);

figure()
plot3(Rglobal(:,1),Rglobal(:,2),Rglobal(:,3))
axis equal
grid on
hold on
plot3(RI(1), RI(2), RI(3),'*m')
plot3(RF(1), RF(2), RF(3),'*b')

%% Non linear interpolator
output2 = NL_interpolator( RI , RF , VI , VF , N_rev , TOF ,M ,PS.Is ,sim );

figure()
subplot(5,1,1)
plot(output2.t*sim.TU/86400,output2.Thrust(:,1));
xlabel('Time [days]')
ylabel('In-plane Thrust [N]')

subplot(5,1,2)
plot(output2.t*sim.TU/86400,180/pi*output2.Thrust(:,2));
xlabel('Time [days]')
ylabel('In-plane Thrust angle [deg]')

subplot(5,1,3)
plot(output2.t*sim.TU/86400,output2.Thrust(:,3));
xlabel('Time [days]')
ylabel('out-of-plane Thrust [N]')

subplot(5,1,4)
plot(output2.t*sim.TU/86400,sqrt(output2.Thrust(:,1).^2 + output2.Thrust(:,3).^2));
xlabel('Time [days]')
ylabel('Thrust [N]')

subplot(5,1,5)
plot(output2.t*sim.TU/86400,output2.m);
xlabel('Time [days]')
ylabel('Mass [kg]')

R_sdrp2  = [output2.r.*cos(output2.theta) output2.r.*sin(output2.theta) output2.z];
Rglobal2 = rotate_local2ecplitic(RI,R_sdrp2,sim.n_sol,output2.Href);

figure()
plot3(Rglobal2(:,1),Rglobal2(:,2),Rglobal2(:,3),'Color',colors(1,:),'DisplayName','Trajectory')
axis equal
grid on
hold on
plot3(RI(1), RI(2), RI(3),'*','Color',colors(2,:),'DisplayName','Start')
plot3(RF(1), RF(2), RF(3),'*','Color',colors(3,:),'DisplayName','End')
legend('show')

%% NLI2
% problem, HI e h_ref sono troppo diversi non so perchè
output3 = NLI2( RI , RF , VI , VF , N_rev , TOF ,M ,PS.Is ,sim );
figure()
subplot(5,1,1)
plot(output3.t*sim.TU/86400,output3.Thrust(:,1));
xlabel('Time [days]')
ylabel('In-plane Thrust [N]')

subplot(5,1,2)
plot(output3.t*sim.TU/86400,180/pi*output3.Thrust(:,2));
xlabel('Time [days]')
ylabel('In-plane Thrust angle [deg]')

subplot(5,1,3)
plot(output3.t*sim.TU/86400,output3.Thrust(:,3));
xlabel('Time [days]')
ylabel('out-of-plane Thrust [N]')

subplot(5,1,4)
plot(output3.t*sim.TU/86400,sqrt(output3.Thrust(:,1).^2 + output3.Thrust(:,3).^2));
xlabel('Time [days]')
ylabel('Thrust [N]')

subplot(5,1,5)
plot(output3.t*sim.TU/86400,output3.m);
xlabel('Time [days]')
ylabel('Mass [kg]')

R_sdrp3  = [output3.r.*cos(output3.theta) output3.r.*sin(output3.theta) output3.z];
Rglobal3 = rotate_local2ecplitic(RI,R_sdrp3,sim.n_sol,output3.Href);

figure()
plot3(Rglobal3(:,1),Rglobal3(:,2),Rglobal3(:,3),'Color',colors(1,:),'DisplayName','Trajectory')
axis equal
grid on
hold on
plot3(RI(1), RI(2), RI(3),'*','Color',colors(2,:),'DisplayName','Start')
plot3(RF(1), RF(2), RF(3),'*','Color',colors(3,:),'DisplayName','End')
legend('show')

TOF
output.t(end)
output2.t(end)
output3.t(end)