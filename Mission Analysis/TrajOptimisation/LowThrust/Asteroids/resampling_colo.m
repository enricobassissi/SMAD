%% RESAMPLING
time_old_1 = sol.departure_mjd2000+output.t1*sim.TU/86400;
dt_new = 0.2; % day, fine discretisation
N_old = length(time_old_1);
dt_old_1 = (time_old_1(end)-time_old_1(1))/N_old; % it's not constant spaced lol
multipl_new = ceil(dt_old_1/dt_new);
N_new = N_old*multipl_new;
time_new_1 = linspace(time_old_1(1),time_old_1(end),N_new);
RTO1_new = interp1(time_old_1,R_transf_orbit_1,time_new_1,'linear');

figure()
plot(time_old_1,R_transf_orbit_1,'DisplayName','Old')
hold on
plot(time_new_1,RTO1_new,'DisplayName','New')
legend('show')

figure()
plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
    'Color',colors(1,:),'DisplayName','Traj SC1')
figure()
plot3(RTO1_new(:,1),RTO1_new(:,2),RTO1_new(:,3),...
    'Color',colors(1,:),'DisplayName','Traj SC1')

%% orbit plots
% transfer orbits
r_transf_orbit_1  = [output.r.leg1.*cos(output.theta.leg1), ...
    output.r.leg1.*sin(output.theta.leg1), output.z.leg1];
R_transf_orbit_1 = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_1,sim.n_sol,output.Href.leg1);


r_transf_orbit_2  = [output.r.leg2.*cos(output.theta.leg2), ...
    output.r.leg2.*sin(output.theta.leg2), output.z.leg2];
R_transf_orbit_2 = rotate_local2ecplitic(r_encounter.astD1,r_transf_orbit_2,sim.n_sol,output.Href.leg2);

% Coasting on Ast 1
start_ast1_sec = (sol.departure_mjd2000+sol.TOF1)*86400;
end_ast1_sec = (sol.departure_mjd2000+sol.TOF1+sol.CT1)*86400;
time_int23 = linspace(start_ast1_sec,end_ast1_sec,50);
y0_coast_ast1 = [r_encounter.astA1*sim.DU; v_encounter.astA1*sim.DU/sim.TU]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[tC1,yC1] = ode113(@rates, time_int23,y0_coast_ast1,options,'sun');

figure()
plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
    'Color',colors(1,:),'DisplayName','Traj SC1')
hold on
hc1 = plot3( yC1(:,1)/sim.DU, yC1(:,2)/sim.DU, yC1(:,3)/sim.DU,'*','Color',colors(1,:),...
    'Markersize',3,'DisplayName','CT');
hc1.Annotation.LegendInformation.IconDisplayStyle = 'off';
hpt2 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
    'Color',colors(2,:));
% hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Earth Dep')
plot3(r_encounter.astA1(1),r_encounter.astA1(2),r_encounter.astA1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1+' Arr')
plot3(r_encounter.astD1(1),r_encounter.astD1(2),r_encounter.astD1(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_1+' Dep')
plot3(r_encounter.astA2(1),r_encounter.astA2(2),r_encounter.astA2(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2+' Arr')

axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
% PLANETS
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
plot_planet_orbit(x(1)*sim.TU/(3600*24),4,colors,6); % mars
% Asteroids
fraction_of_the_orbit = 1;
hello_orbit1 = sol.departure_mjd2000+output.t1(end)*sim.TU/(3600*24);
% hello_orbit2 = sol.departure_mjd2000+(output.t1(end)+output.t2(end))*sim.TU/(3600*24);
plot_asteorid_orbit(hello_orbit1,fraction_of_the_orbit,sol.asteroid_1,colors,3);
% plot_asteorid_orbit(hello_orbit2,fraction_of_the_orbit,sol.asteroid_2,colors,4);

plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)
