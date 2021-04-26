%% simulation parameters
sim.mu    = 132712440018          ; % actractor parameter [km^3 s^-2]
sim.DU    = 149597870.7           ; % distance unit [km]
sim.TU    = (sim.DU^3/sim.mu )^0.5; % time unit [s]
sim.mu    = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol = 100;                    % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 
sim.out_shape = 2;                  % out-of-plane shape
sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = -1;                     % direction of integration (1 FW, 2 BW)

%% Propulsive system parameters
PS.Is = 3000/sim.TU;  % non-dimensional specific impulse


%% 
RI = [ 1  0 0.1]'; % initial position [DU]
RF = [-1 -1 0]'; % final position [DU]
VI = [ 0  1 0]'; % initial velocity [DU/TU]
VF = [ 0.7  -0.7 0]'; % final velocity [DU/TU]
N_rev = 2; % number of revolution
TOF = 25.4; % TOF [TU]
M = 1000; % SC mass [kg]
hp = 3; 
kp = 3;

[ output] = CW_LowLambert( RI , RF , VI , VF , N_rev , TOF ,M ,hp , kp , PS ,sim )

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
plot3(output.r.*cos(output.l),output.r.*sin(output.l),output.z)
axis equal
grid on