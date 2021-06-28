% --- alessio stuff
data.n_int = 50;
data.MU = sol.M_start_SC1_leg1; % mass adimensionalisation on wet mass
SC1.asteroid_1 = sol.asteroid_1;

v_launcher = sol.v_inf_magn/sim.DU*sim.TU*[cos(sol.el)*cos(sol.az); cos(sol.el)*sin(sol.az); sin(sol.el)];
v_dep = v_encounter.EA + v_launcher;  %if parabolic escape (v_extra = 0)
SC1.leg1.v_in = v_dep;
SC1.leg1.v_end = v_encounter.astA1;
SC1.leg1.r_in = r_encounter.EA;
SC1.leg1.r_end = r_encounter.astA1;
SC1.leg1.N_rev = sol.Nrev(1);
SC1.leg1.TOF = sol.TOF1_ADIM;

DT_input = SC1.leg1;
[SC1.leg1.HS, SC1.leg1.ANGLES_AS, SC1.leg1.timead, SC1.leg1.Xad, ...
    SC1.leg1.Xpropad_DT, SC1.leg1.Xpropad_AS] = DT_executable_attitude_sensitivity(DT_input, sim, data);

% ---- cartesian orbit
% ---- nominal propagation
r_transf_orbit_nom  = [SC1.leg1.HS.X(:,1).*cos(SC1.leg1.HS.X(:,2)), ...
    SC1.leg1.HS.X(:,1).*sin(SC1.leg1.HS.X(:,2)), SC1.leg1.HS.X(:,3)];
R_cartesian_nom = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_nom,data.n_int,sol.Href(:,1));
% ---- DT propagation
r_transf_orbit_DT  = [SC1.leg1.Xpropad_DT(:,1).*cos(SC1.leg1.Xpropad_DT(:,2)), ...
    SC1.leg1.Xpropad_DT(:,1).*sin(SC1.leg1.Xpropad_DT(:,2)), SC1.leg1.Xpropad_DT(:,3)];
R_cartesian_DT = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_DT,data.n_int,sol.Href(:,1));
% ---- AS propagation
r_transf_orbit_AS  = [SC1.leg1.Xpropad_AS(:,1).*cos(SC1.leg1.Xpropad_AS(:,2)), ...
    SC1.leg1.Xpropad_AS(:,1).*sin(SC1.leg1.Xpropad_AS(:,2)), SC1.leg1.Xpropad_AS(:,3)];
R_cartesian_AS = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_AS,data.n_int,sol.Href(:,1));

% ----- plots
figure('Name','X comparison')
% plot(SC1.leg1.timead, SC1.leg1.Xad,'Color',colors(1,:),'DisplayName','Xad');
hold on; grid on;
hp1 = plot(SC1.leg1.timead*sim.TU/86400, SC1.leg1.Xpropad_DT,'Color',colors(2,:));
hp2 = plot(SC1.leg1.timead*sim.TU/86400, SC1.leg1.Xpropad_AS,'Color',colors(3,:));
xlabel(' time [d]'); ylabel('dynamics'); % legend([hp1, hp2], 'DT Prop','DT Prop Att');

figure('Name','angles')
plot(SC1.leg1.timead*sim.TU/86400, rad2deg(SC1.leg1.HS.beta),'Color',colors(1,:),'DisplayName','HS el');
hold on; grid on;
plot(SC1.leg1.timead*sim.TU/86400, rad2deg(SC1.leg1.HS.alpha),'Color',colors(1,:),'DisplayName','HS az');
plot(SC1.leg1.timead*sim.TU/86400, rad2deg(SC1.leg1.ANGLES_AS(:,1)),'Color',colors(2,:),'DisplayName','AS az');
plot(SC1.leg1.timead*sim.TU/86400, rad2deg(SC1.leg1.ANGLES_AS(:,2)),'Color',colors(2,:),'DisplayName','HS az');
legend('show'); xlabel(' time [d]'); ylabel('angles [deg]');

figure('Name','Angles rand difference')
plot(SC1.leg1.timead*sim.TU/86400, rad2deg(SC1.leg1.HS.beta - SC1.leg1.ANGLES_AS(:,2)),'Color',colors(1,:),'DisplayName','el diff');
hold on; grid on;
plot(SC1.leg1.timead*sim.TU/86400, rad2deg(SC1.leg1.HS.alpha - SC1.leg1.ANGLES_AS(:,1)),'Color',colors(2,:),'DisplayName','az diff');
legend('show'); xlabel(' time [d]'); ylabel('rand: nominal - AS [deg]');

figure('Name','max err on traj cartesian coordinates')
plot(SC1.leg1.timead*sim.TU/86400, (R_cartesian_DT(:,1) - R_cartesian_AS(:,1))*sim.DU,'Color',colors(1,:),'DisplayName','x diff');
hold on; grid on;
plot(SC1.leg1.timead*sim.TU/86400, (R_cartesian_DT(:,2) - R_cartesian_AS(:,2))*sim.DU,'Color',colors(2,:),'DisplayName','y diff');
plot(SC1.leg1.timead*sim.TU/86400, (R_cartesian_DT(:,3) - R_cartesian_AS(:,3))*sim.DU,'Color',colors(3,:),'DisplayName','z diff');
plot(SC1.leg1.timead*sim.TU/86400, (vecnorm(R_cartesian_DT,2,2) - vecnorm(R_cartesian_AS,2,2))*sim.DU,'Color',colors(4,:),'DisplayName','norm diff');
legend('show'); xlabel(' time [d]'); ylabel('positions [km]');

figure('Name','3D plot orbit')
plot3(R_cartesian_DT(:,1),R_cartesian_DT(:,2),R_cartesian_DT(:,3),'Color',colors(1,:),'DisplayName','DT orbit');
hold on; grid on; axis equal;
plot3(R_cartesian_DT(end,1),R_cartesian_DT(end,2),R_cartesian_DT(end,3),'o','Color',colors(1,:),'DisplayName','DT end');
plot3(R_cartesian_AS(:,1),R_cartesian_AS(:,2),R_cartesian_AS(:,3),'--','Linewidth',3.5,'Color',colors(2,:),'DisplayName','AS orbit');
plot3(R_cartesian_AS(end,1),R_cartesian_AS(end,2),R_cartesian_AS(end,3),'o','Color',colors(2,:),'DisplayName','AS end');
legend('show'); xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');