
clear; close all; clc;

%% Adimensional position and velocity 
r1 = [0.717953940614415;-0.716259972545788;0];
v1 = [0.690130064370903;0.704151274121327;0];

r2 =  [-1.560337259616491;-0.460895252843023;0.028886386394797];
v2 = [0.260165078626092;-0.710377652565987;-0.021219637204627];

%% Simulation parameters
sim.mu_dim    = 132712440018          ; % actractor parameter [km^3 s^-2]
sim.DU        = 149597870.7           ; % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu        = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol     = 200;                    % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 

sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = 1;                     % direction of integration (1 FW, -1 BW)


sim.PS.Isp = 3000/sim.TU;  % non-dimensional specific impulse

sim.M = 1000; % SC mass [kg]

TOF1 = 25;
sim.TOF_imposed_flag = 0; % non impone TOF

[output] = NL_interpolator( r1 , r2 , v1 , v2 , 0 , TOF1 , sim.M ,sim.PS.Isp ,sim );

figure()
subplot(4,1,1)
plot(output.t*sim.TU/86400,output.Thrust(:,1));
xlabel('Time [days]')
ylabel('In-plane Thrust [N]')

subplot(4,1,2)
plot(output.t*sim.TU/86400,180/pi*output.Thrust(:,2));
xlabel('Time [days]')
ylabel('In-plane Thrust angle [deg]')

subplot(4,1,3)
plot(output.t*sim.TU/86400,output.Thrust(:,3));
xlabel('Time [days]')
ylabel('out-of-plane Thrust [N]')


figure;
plot(output.t*sim.TU/86400,sqrt(output.Thrust(:,1).^2 + output.Thrust(:,3).^2));
xlabel('Time [days]')
ylabel('Thrust [N]')

figure;
plot(output.t*sim.TU/86400,output.m);
xlabel('Time [days]')
ylabel('Mass [kg]')

r3  = [output.r.*cos(output.theta) output.r.*sin(output.theta) output.z];
R3 = rotate_local2ecplitic(r1,r3,sim.n_sol,output.Href);

day1 = [2028 1 1 0 0 0];
day2 = [2031 1 1 0 0 0];

t1 = date2mjd2000(day1);
t2 = date2mjd2000(day2);
times = linspace(t1,t2,1000);

for i=1:length(times)
    % Orbit 1
    [kep1,ksun] = uplanet(times(i),3);
    [ra(i,1:3),va(i,1:3)] = sv_from_coe(kep1,ksun);
    ra(i,1:3) = ra(i,1:3)/sim.DU;
    
    % Orbit 2
    [kep2,ksun] = uplanet(times(i),4);
    [rb(i,1:3),~] = sv_from_coe(kep2,ksun);
    rb(i,1:3) = rb(i,1:3)/sim.DU;
end

figure()
plot3(R3(:,1),R3(:,2),R3(:,3))
hold on
plot3(r1(1),r1(2),r1(3),'*m')
plot3(r2(1),r2(2),r2(3),'*c')
axis equal
grid on
hold on
% Earth
hold on
plot3(ra(:,1),ra(:,2),ra(:,3),'--g'); % geocentric equatorial frame 
% Mars
hold on
plot3(rb(:,1),rb(:,2),rb(:,3),'--r');
% Sun
plot3(0,0,0,'oy')
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
legend('show')

