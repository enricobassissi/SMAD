%% ---------- transfer orbits Positions --------------- %%
%% SC1
% EA Dep -> DSM1
start_EA_DSM1_sec = (sol.departure_mjd2000)*86400;
end_EA_DSM_sec_SC1 = (sol.departure_mjd2000+sol.TOF_DSM1)*86400;
time_EA_DSM_SC1 = linspace(start_EA_DSM1_sec,end_EA_DSM_sec_SC1,100);
y0_EA_DSM_SC1 = [r_encounter.EA*sim.DU; v_encounter.EA_DSM1*sim.DU/sim.TU]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_EA_DSM_SC1,y_EA_DSM_SC1] = ode113(@rates, time_EA_DSM_SC1,y0_EA_DSM_SC1,options,'sun');

% DSM1 -> GA1
start_DSM_GA_sec_SC1 = (sol.departure_mjd2000+sol.TOF_DSM1)*86400;
end_DSM_GA_sec_SC1 = (sol.departure_mjd2000+sol.TOF_DSM1+sol.TOFGA_1)*86400;
time_DSM_GA_SC1 = linspace(start_DSM_GA_sec_SC1,end_DSM_GA_sec_SC1,100);
y0_DSM_GA_SC1 = [r_encounter.DSM1*sim.DU; v_encounter.DSM1_GA1*sim.DU/sim.TU]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_DSM_GA_SC1,y_DSM_GA_SC1] = ode113(@rates, time_DSM_GA_SC1,y0_DSM_GA_SC1,options,'sun');

% GA1 -> Ast 1 arr
r_transf_orbit_1  = [output.r.leg1.*cos(output.theta.leg1), ...
    output.r.leg1.*sin(output.theta.leg1), output.z.leg1];
R_transf_orbit_1 = rotate_local2ecplitic(r_encounter.GA1,r_transf_orbit_1,sim.n_sol,output.Href.leg1);

% Ast 1 dep -> Ast 2 arr
r_transf_orbit_2  = [output.r.leg2.*cos(output.theta.leg2), ...
    output.r.leg2.*sin(output.theta.leg2), output.z.leg2];
R_transf_orbit_2 = rotate_local2ecplitic(r_encounter.astD1,r_transf_orbit_2,sim.n_sol,output.Href.leg2);

%% SC2
% EA Dep -> DSM2
start_EA_DSM2_sec = (sol.departure_mjd2000)*86400;
end_EA_DSM2_sec_SC2 = (sol.departure_mjd2000+sol.TOF_DSM2)*86400;
time_EA_DSM2_SC2 = linspace(start_EA_DSM2_sec,end_EA_DSM2_sec_SC2,100);
y0_EA_DSM2_SC2 = [r_encounter.EA*sim.DU; v_encounter.EA_DSM2*sim.DU/sim.TU]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_EA_DSM2_SC2,y_EA_DSM2_SC2] = ode113(@rates, time_EA_DSM2_SC2,y0_EA_DSM2_SC2,options,'sun');

% DSM2 -> GAa
start_DSM2_GAa_sec_SC2 = (sol.departure_mjd2000+sol.TOF_DSM2)*86400;
end_DSM2_GAa_sec_SC2 = (sol.departure_mjd2000+sol.TOF_DSM2+sol.TOFGA_a)*86400;
time_DSM2_GAa_SC2 = linspace(start_DSM2_GAa_sec_SC2,end_DSM2_GAa_sec_SC2,100);
y0_DSM2_GAa_SC2 = [r_encounter.DSM2*sim.DU; v_encounter.DSM2_GAa*sim.DU/sim.TU]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_DSM2_GAa_SC2,y_DSM2_GAa_SC2] = ode113(@rates, time_DSM2_GAa_SC2,y0_DSM2_GAa_SC2,options,'sun');

% GAa -> Ast a arr
r_transf_orbit_a  = [output.r.lega.*cos(output.theta.lega), ...
    output.r.lega.*sin(output.theta.lega), output.z.lega];
R_transf_orbit_a = rotate_local2ecplitic(r_encounter.GAa,r_transf_orbit_a,sim.n_sol,output.Href.lega);

% % Coasting on Ast a
% start_asta_sec = (sol.departure_mjd2000+sol.TOFa)*86400;
% end_asta_sec = (sol.departure_mjd2000+sol.TOFa+sol.CTa)*86400;
% % time_int_SC2 = linspace(start_asta_sec,end_asta_sec,50);
% time_int_SC2 = sol.departure_mjd2000*86400+output.CTa*sim.TU;
% y0_coast_asta = [r_encounter.astAa*sim.DU; v_encounter.astAa*sim.DU/sim.TU]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
% options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
% [tCa,yCa] = ode113(@rates, time_int_SC2, y0_coast_asta,options,'sun');

% Ast a dep -> Ast b arr
r_transf_orbit_b  = [output.r.legb.*cos(output.theta.legb), ...
    output.r.legb.*sin(output.theta.legb), output.z.legb];
R_transf_orbit_b = rotate_local2ecplitic(r_encounter.astDa,r_transf_orbit_b,sim.n_sol,output.Href.legb);


%% Plot SC1
figure('Name','SC 1 Travel')
hl1 = plot3( y_EA_DSM_SC1(:,1)/sim.DU, y_EA_DSM_SC1(:,2)/sim.DU, y_EA_DSM_SC1(:,3)/sim.DU,...
    'Color',colors(1,:),'DisplayName','EA - DSM');
hold on
hl2 = plot3( y_DSM_GA_SC1(:,1)/sim.DU, y_DSM_GA_SC1(:,2)/sim.DU, y_DSM_GA_SC1(:,3)/sim.DU,...
    'Color',colors(2,:),'DisplayName','DSM - GA');
hlt1 = plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
    'Color',colors(3,:),'DisplayName','GA - Ast1');
% hlt1.Annotation.LegendInformation.IconDisplayStyle = 'off';
hlt2 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
    'Color',colors(4,:),'DisplayName','Ast1 - Ast2');
% hlt2.Annotation.LegendInformation.IconDisplayStyle = 'off';

% ---- points ------- %
plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Earth Dep')
plot3(r_encounter.DSM1(1),r_encounter.DSM1(2),r_encounter.DSM1(3),...
    'x','Color',colors(1,:),'DisplayName','DSM')
plot3(r_encounter.GA1(1),r_encounter.GA1(2),r_encounter.GA1(3),...
    'd','Color',colors(8,:),'DisplayName','Earth GA')
plot3(r_encounter.astA1(1),r_encounter.astA1(2),r_encounter.astA1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1+' Arr')
plot3(r_encounter.astD1(1),r_encounter.astD1(2),r_encounter.astD1(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_1+' Dep')
plot3(r_encounter.astA2(1),r_encounter.astA2(2),r_encounter.astA2(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2+' Arr')

% ----- PLANETS and SUN ---- %
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
plot_planet_orbit(x(1)*sim.TU/(3600*24),4,colors,6); % mars
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')

% ---- Asteroids ---- %
fraction_of_the_orbit = 1;
hello_orbit1 = sol.departure_mjd2000+sol.TOF1;
hello_orbit2 = sol.departure_mjd2000+(sol.TOF1+sol.CT1+sol.TOF2);
plot_asteorid_orbit(hello_orbit1,fraction_of_the_orbit,sol.asteroid_1,colors,3);
plot_asteorid_orbit(hello_orbit2,fraction_of_the_orbit,sol.asteroid_2,colors,4);

% ---- options ----- %
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
clearvars hl1 hl2 hlt1 hlt2
legend('show')
view(2)

%% Plots SC2
figure('Name','SC 2 Travel')
hl1 = plot3( y_EA_DSM2_SC2(:,1)/sim.DU, y_EA_DSM2_SC2(:,2)/sim.DU, y_EA_DSM2_SC2(:,3)/sim.DU,...
    'Color',colors(1,:),'DisplayName','EA - DSM');
hold on
hl2 = plot3( y_DSM2_GAa_SC2(:,1)/sim.DU, y_DSM2_GAa_SC2(:,2)/sim.DU, y_DSM2_GAa_SC2(:,3)/sim.DU,...
    'Color',colors(2,:),'DisplayName','DSM - GA');
hlt1 = plot3(R_transf_orbit_a(:,1),R_transf_orbit_a(:,2),R_transf_orbit_a(:,3),...
    'Color',colors(3,:),'DisplayName','GA - Ast a');
% % hca = plot3( yCa(:,1)/sim.DU, yCa(:,2)/sim.DU, yCa(:,3)/sim.DU,'*','Color',colors(2,:),...
% %     'Markersize',3);
% % hca.Annotation.LegendInformation.IconDisplayStyle = 'off';
hlt2 = plot3(R_transf_orbit_b(:,1),R_transf_orbit_b(:,2),R_transf_orbit_b(:,3),...
    'Color',colors(4,:),'DisplayName','Ast a - Ast b');
% hlt2.Annotation.LegendInformation.IconDisplayStyle = 'off';

% ---- points ----- %
plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Earth Dep')
plot3(r_encounter.DSM2(1),r_encounter.DSM2(2),r_encounter.DSM2(3),...
    'x','Color',colors(1,:),'DisplayName','DSM')
plot3(r_encounter.GAa(1),r_encounter.GAa(2),r_encounter.GAa(3),...
    'd','Color',colors(8,:),'DisplayName','Earth GA')
plot3(r_encounter.astAa(1),r_encounter.astAa(2),r_encounter.astAa(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_a+' Arr')
plot3(r_encounter.astDa(1),r_encounter.astDa(2),r_encounter.astDa(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_a+' Dep')
plot3(r_encounter.astAb(1),r_encounter.astAb(2),r_encounter.astAb(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_b+' Arr')
% ----- PLANETS and SUN ---- %
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
plot_planet_orbit(x(1)*sim.TU/(3600*24),4,colors,6); % mars
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')

% ---- Asteroids ---- %
fraction_of_the_orbit = 1;
hello_orbita = sol.departure_mjd2000+sol.TOFa;
hello_orbitb = sol.departure_mjd2000+sol.TOFa+sol.CTa+sol.TOFb;
plot_asteorid_orbit(hello_orbita,fraction_of_the_orbit,sol.asteroid_a,colors,3);
plot_asteorid_orbit(hello_orbitb,fraction_of_the_orbit,sol.asteroid_b,colors,4);

% ---- options ----- %
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
clearvars hl1 hl2 hlt1 hlt2
legend('show')
view(2)