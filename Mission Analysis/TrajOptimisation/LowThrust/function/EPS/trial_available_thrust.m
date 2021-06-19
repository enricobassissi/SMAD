addpath(genpath(pwd))
clc,clear all,close all

load('trajectories_2SC.mat');
load('power_propulsion_data.mat');
[~,~,~,~,~,TT] = saa_estimation(Traj_SC1);
t = linspace(0,5,length(Traj_SC1))';

% power_propulsion_data.P_other_subsystem_margined= 100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to overwrite

T_available = available_thrust(t, Traj_SC1, TT, power_propulsion_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%distance from the sun
d_sun = vecnorm(Traj_SC1,2,2);
eta_space = solar_irradiance_cooler(d_sun, 1, 1);

% angle of panels respect to sun
a_sun = acos(dot(Traj_SC1,TT,2)./(vecnorm(TT,2,2).*d_sun)) -pi/2; 
eta_angle = cos(a_sun);

% degradation
eta_degradation = (1-power_propulsion_data.D_cell_radiation-power_propulsion_data.D_cell_other).^t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, hold on
yyaxis left
plot(t,eta_space)
plot(t,eta_angle)
plot(t,eta_degradation)
grid on
ylabel('\eta_x')
yyaxis right
plot(t,T_available)
xlabel('Time [y]')
ylabel('T_{available} [mN]')
legend('distance','angle','degradation','thrust','Location','Best')