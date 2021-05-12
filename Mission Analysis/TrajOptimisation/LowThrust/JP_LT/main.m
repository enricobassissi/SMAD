

clc;clear;close all;

%% simulation parameters
sim.mu_dim    = 132712440018          ; % actractor parameter [km^3 s^-2] -- sun
sim.DU    = 149597870.7           ; % distance unit [km]
sim.TU    = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu    = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol = 100;                    % number of computational nodes %% 100 initially
sim.x = linspace(0,1,sim.n_sol)';   % 
sim.out_shape = 2;                  % out-of-plane shape (2:CONWAY)
sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = -1;                     % direction of integration (1 FW, -1 BW)

%% Propulsive system parameters
PS.Is = 3000/sim.TU;  % non-dimensional specific impulse


%% 
RI = [ 1  0 0.1]'; % initial position [DU]  %% ORIGINALE: [1 0 0.1] 
RF = [-1 -1 0]'; % final position [DU]
% RI = [0.9545
%     0.2938
%          0]; % initial position [DU]  %% ORIGINALE: [1 0 0.1] 
% RF = [    1.3767
%    -0.1449
%    -0.0375]; % final position [DU]
VI = [ 0  1 0]'; % initial velocity [DU/TU]
VF = [ 0.7  -0.7 0]'; % final velocity [DU/TU]
% VI = [   -0.3105
%     0.9520
%          0]; % initial velocity [DU/TU]
% VF = [0.1122
%     0.8804
%     0.0151]; % final velocity [DU/TU]

N_rev = 2; % number of revolution
TOF = 25.4; % TOF [TU] %25.4
M = 1000; % SC mass [kg]
hp = 3; %3 OUT OF PLANE SHAPE PARAMETERS
kp = 3; %3

% %%
% [output] = CW_LowLambert( RI , RF , VI , VF , N_rev , TOF ,M ,hp , kp , PS ,sim );
% 
% output.m(1);
% 
% %%
% figure()
% subplot(5,1,1)
% plot(output.t*sim.TU/86400,output.u(:,1));
% xlabel('Time [days]')
% ylabel('In-plane Thrust [N]')
% 
% subplot(5,1,2)
% plot(output.t*sim.TU/86400,180/pi*output.u(:,2));
% xlabel('Time [days]')
% ylabel('In-plane Thrust angle [deg]')
% 
% subplot(5,1,3)
% plot(output.t*sim.TU/86400,output.u(:,3));
% xlabel('Time [days]')
% ylabel('out-of-plane Thrust [N]')
% 
% subplot(5,1,4)
% plot(output.t*sim.TU/86400,sqrt(output.u(:,1).^2 + output.u(:,3).^2));
% xlabel('Time [days]')
% ylabel('Thrust [N]')
% 
% subplot(5,1,5)
% plot(output.t*sim.TU/86400,output.m);
% xlabel('Time [days]')
% ylabel('Mass [kg]')
% 
% R_sdrp  = [output.r.*cos(output.l) output.r.*sin(output.l) output.z];
% Rglobal = rotate_local2ecplitic(RI,R_sdrp,sim.n_sol,output.h_ref_v);
% 
% figure()
% plot3(Rglobal(:,1),Rglobal(:,2),Rglobal(:,3))
% axis equal
% grid on
% hold on
% plot3(RI(1), RI(2), RI(3),'*m')
% plot3(RF(1), RF(2), RF(3),'*b')

%%

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

R_sdrp2  = [output2.r.*cos(output2.theta) output2.r.*sin(output2.theta) output2.z];
Rglobal2 = rotate_local2ecplitic(RI,R_sdrp2,sim.n_sol,output2.Href);

figure()
plot3(Rglobal2(:,1),Rglobal2(:,2),Rglobal2(:,3))
axis equal
grid on
hold on
plot3(RI(1), RI(2), RI(3),'*m')
plot3(RF(1), RF(2), RF(3),'*b')

output.t(end)
output2.t(end)