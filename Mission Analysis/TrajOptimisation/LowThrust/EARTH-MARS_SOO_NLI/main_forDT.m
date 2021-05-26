
clear; close all; clc;

%% Initial and final values
%-- Initial point
r1 = 1;
theta1 = 0;

x1 = r1*cos(theta1); y1 = r1*sin(theta1); z1 = 0;
R1 = [x1; y1; z1];

vr1 = 0; 
vth1 = 1; 

r1_t = vr1;
theta1_t = vth1/r1;

x1_t = r1_t*cos(theta1) - r1*sin(theta1)*theta1_t;
y1_t = r1_t*sin(theta1) + r1*cos(theta1)*theta1_t;
z1_t = 0;

V1 = [x1_t; y1_t; z1_t];

%-- Final point
r2 = 4;
theta2 = 0; %%%

x2 = r2*cos(theta2); y2 = r2*sin(theta2); z2 = 0;
R2 = [x2; y2; z2];

vr2 = 0;
vth2 = 0.5;

r2_t = vr2;
theta2_t = vth2/r2;

x2_t = r2_t*cos(theta2) - r2*sin(theta2)*theta2_t;
y2_t = r2_t*sin(theta2) + r2*cos(theta2)*theta2_t;
z2_t = 0;

V2 = [x2_t; y2_t; z2_t];


%% simulation parameters
sim.mu_dim    = 398600;                     % actractor parameter [km^3 s^-2]
sim.DU        = 6378 ;                      % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu        = 1;                          % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol     = 500;                        % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   

sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = 1;                      % direction of integration (1 FW, -1 BW), 
                                        %  1 means imposing wet mass at beginning

sim.PS.Isp = 3000/sim.TU;  % non-dimensional specific impulse
sim.M = 1000;

sim.TOF_imposed_flag = 0; % 0 means TOF free
TOF = 25; % puoi mettere un valore a caso, tanto siamo nel caso TOF FREE

N = 6;

sim.out_shape = 2; %out-of-plane shape (2:CONWAY)
sim.tol_vers = 1e-8;

PS.Is = 3000/sim.TU;
% [output] = NL_interpolator( R1 , R2 , V1 , V2 , N , TOF , sim.M ,sim.PS.Isp ,sim );
% 
% r3  = [output.r.*cos(output.theta) output.r.*sin(output.theta) output.z];
% R3 = rotate_local2ecplitic(R1,r3,sim.n_sol,output.Href);

[ output ] = CW_LowLambert( R1 , R2 , V1 , V2 , N , TOF ,sim.M ,3 , 3 , PS ,sim )
r3  = [output.r.*cos(output.l) output.r.*sin(output.l) output.z];
R3 = rotate_local2ecplitic(R1,r3,sim.n_sol,output.h_ref_v);

figure()
plot3(R3(:,1),R3(:,2),R3(:,3))
hold on
plot3(R1(1),R1(2),R1(3),'*m')
plot3(R2(1),R2(2),R2(3),'*c')
axis equal
grid on
hold on
xlabel('x [DU]'); ylabel('y [DU]'); ylabel('y [DU]'); 
legend('show')

figure;
plot(output.t*sim.TU/86400,sqrt(output.u(:,1).^2 + output.u(:,3).^2));
xlabel('Time [days]')
ylabel('Thrust Dimensional [N]')

figure;
plot(output.t*sim.TU/86400,sqrt(output.u(:,1).^2 + output.u(:,3).^2)./output.m*sim.TU^2/sim.DU);
xlabel('Time [days]')
ylabel('Acceleration Adimensional [-]')

amax = 0.01;

figure;
plot(output.t*sim.TU/86400,output.m);
xlabel('Time [days]')
ylabel('Mass [kg]')
