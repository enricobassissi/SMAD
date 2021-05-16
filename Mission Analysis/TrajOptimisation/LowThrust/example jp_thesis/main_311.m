clc;clear;close all;

%% simulation parameters
sim.mu_dim  = astroConstants(20)        ; % actractor parameter [km^3 s^-2] -- sun
sim.DU      = 1738;                       % distance unit [km]
sim.TU      = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu      = 1;                          % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol   = 100;                        % number of computational nodes %% 100 initially
sim.x = linspace(0,1,sim.n_sol)';         % 
sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU));   % non-dimensional g0
sim.direction = -1;                       % direction of integration (1 FW, -1 BW)

%% Propulsive system parameters
PS.Is = 3800/sim.TU;  % non-dimensional specific impulse


%% 
coe1 = [4.5*sim.DU 0.2 deg2rad(5) 0 deg2rad(20) deg2rad(10)];
coe2 = [7*sim.DU 0.1 deg2rad(15) deg2rad(10) deg2rad(5) deg2rad(350)];

T1 = 2*pi*sqrt((4.5*sim.DU)^3/astroConstants(20));
T2 = 2*pi*sqrt((7*sim.DU)^3/astroConstants(20));


[RI VI] = sv_from_coe(coe1,astroConstants(20));
RI = RI/sim.DU;
VI = VI/sim.DU*sim.TU;

[RF VF] = sv_from_coe(coe2,astroConstants(20));
RF = RF/sim.DU;
VF = VF/sim.DU*sim.TU;


N_rev = 0; % number of revolution
TOF = 24.02*3600/sim.TU; % TOF [TU] 
M = 1000; % SC mass [kg]

IC1 = [RI*sim.DU; VI*sim.DU/sim.TU];
IC2 = [RF*sim.DU; VF*sim.DU/sim.TU];

options= odeset('RelTol', 1e-13, 'AbsTol',1e-14);
[T11 R11] = ode113(@rates,[0 T1] ,IC1,options,'moon');

[T22 R22] = ode113(@rates,[0 T2] ,IC2,options,'moon');

sim.out_shape = 2;                  % out-of-plane shape (2:CONWAY)
hp = 3; %3 OUT OF PLANE SHAPE PARAMETERS
kp = 3; %3

%% CONWAY
[output] = CW_LowLambert( RI , RF , VI , VF , N_rev , TOF ,M ,hp , kp , PS,sim );

output.m(1)

%%
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
Rglobal = rotate_local2ecplitic(RI,R_sdrp,sim.n_sol,output.h_ref_v)

figure()
plot3(Rglobal(:,1),Rglobal(:,2),Rglobal(:,3))
axis equal
grid on
hold on
plot3(RI(1), RI(2), RI(3),'*m')
plot3(RF(1), RF(2), RF(3),'*b')

fuel_mass1 = output.m(1)-output.m(end)
max_thrust1 = max(sqrt(output.u(:,1).^2 + output.u(:,3).^2))
TOF_hours1 = output.t(end)*sim.TU/3600


%% NLI
output2 = NL_interpolator( RI , RF , VI , VF , N_rev , TOF ,M ,PS.Is ,sim );
 

figure()
subplot(2,3,1)
plot(rad2deg(output2.theta),output2.Thrust(:,1));
grid on
xlabel('Theta[deg]')
ylabel('In-plane Thrust [N]')

subplot(2,3,2)
plot(rad2deg(output2.theta),180/pi*output2.Thrust(:,2));
grid on
xlabel('Theta [deg]')
ylabel('In-plane Thrust angle [deg]')

subplot(2,3,3)
plot(rad2deg(output2.theta),output2.Thrust(:,3));
grid on
xlabel('Theta [deg]')
ylabel('out-of-plane Thrust [N]')

subplot(2,3,4)
plot(rad2deg(output2.theta),sqrt(output2.Thrust(:,1).^2 + output2.Thrust(:,3).^2));
grid on
xlabel('Theta [deg]')
ylabel('Thrust [N]')

subplot(2,3,5)
plot(rad2deg(output2.theta),output2.m);
grid on
xlabel('Theta [deg]')
ylabel('Mass [kg]')

R_sdrp2  = [output2.r.*cos(output2.theta) output2.r.*sin(output2.theta) output2.z];
Rglobal2 = rotate_local2ecplitic(RI,R_sdrp2,sim.n_sol,output2.Href);

figure()
plot3(Rglobal2(:,1),Rglobal2(:,2),Rglobal2(:,3))
axis equal
grid on
hold on
plot3(RI(1), RI(2), RI(3),'*r')
plot3(RF(1), RF(2), RF(3),'*b')
plot3(R11(:,1)/sim.DU,R11(:,2)/sim.DU,R11(:,3)/sim.DU,'--m')
plot3(R22(:,1)/sim.DU,R22(:,2)/sim.DU,R22(:,3)/sim.DU,'--c')


fuel_mass = output2.m(1)-output2.m(end)
max_thrust = max(sqrt(output2.Thrust(:,1).^2 + output2.Thrust(:,3).^2))
TOF_hours = output2.t(end)*sim.TU/3600