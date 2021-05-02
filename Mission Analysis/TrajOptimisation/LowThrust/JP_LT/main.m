

clc;clear;close all;

%% simulation parameters
sim.mu    = 132712440018          ; % actractor parameter [km^3 s^-2] -- sun
sim.DU    = 149597870.7           ; % distance unit [km]
sim.TU    = (sim.DU^3/sim.mu )^0.5; % time unit [s]
sim.mu    = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol = 100;                    % number of computational nodes %% 100 initially
sim.x = linspace(0,1,sim.n_sol)';   % 
sim.out_shape = 2;                  % out-of-plane shape (2:CONWAY)
sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = -1;                     % direction of integration (1 FW, 2 BW)

%% Propulsive system parameters
PS.Is = 3000/sim.TU;  % non-dimensional specific impulse


%% 
RI = [ 1  0 0.1]'; % initial position [DU]  %% ORIGINALE: [1 0 0.1] 
RF = [-1 -1 0]'; % final position [DU]
VI = [ 0  1 0]'; % initial velocity [DU/TU]
VF = [ 0.7  -0.7 0]'; % final velocity [DU/TU]
N_rev = 2; % number of revolution
TOF = 25.4; % TOF [TU]
M = 1000; % SC mass [kg]
hp = 3; %3 OUT OF PLANE SHAPE PARAMETERS
kp = 3; %3

[output] = CW_LowLambert( RI , RF , VI , VF , N_rev , TOF ,M ,hp , kp , PS ,sim );

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

figure()
plot3(output.r.*cos(output.l),output.r.*sin(output.l), output.z)
axis equal
grid on
hold on
plot3(RI(1), RI(2), RI(3),'*m')
plot3(RF(1), RF(2), RF(3),'*b')

%%
% ksun =  astroConstants(4);
% 
% MJD01 = date2mjd2000([2021 1 1 0 0 0]);
% [kep_EA,ksun] = uplanet(MJD01, 3);
% [RI, VI] = sv_from_coe(kep_EA,ksun);
% 
% RI = RI/sim.DU;
% VI = VI/sim.DU*sim.TU;
% 
% % arrival on mars
% MJDF1 = date2mjd2000([2023 1 1 0 0 0]);
% [kep_MA,ksun] = uplanet(MJDF1, 4);
% [RF, VF] = sv_from_coe(kep_MA,ksun);
% 
% RF = RF/sim.DU;
% VF = VF/sim.DU*sim.TU;

Isp = 3000/sim.TU; 
output2 = NL_interpolator( RI , RF , VI , VF , N_rev , TOF ,M ,Isp ,sim );
 
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

figure()
plot3(output2.r.*cos(output2.theta),output2.r.*sin(output2.theta), output2.z)
axis equal
grid on
hold on
plot3(RI(1), RI(2), RI(3),'*m')
plot3(RF(1), RF(2), RF(3),'*b')
