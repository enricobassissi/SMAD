function [sol] = plot_impulsive_soo_DSM(sol,asteroid_names_sequence,data,sim,colors)
    
    AU = astroConstants(2);

    %% initialise stuff from solution
    MJD01 = sol.MJD0;
    MJDF1 = MJD01 + sol.TOF1a + sol.TOF1b;
    MJD02 = MJDF1 + sol.buffer_time1;
    MJDF2 = MJD02 + sol.TOF2;
    MJD03 = MJDF2 + sol.buffer_time2;
    MJDF3 = MJD03 + sol.TOF3;
    MJD04 = MJDF3 + sol.buffer_time3;
    MJDF4 = MJD04 + sol.TOF4;
    
    ast1 = asteroid_names_sequence(1);
    ast2 = asteroid_names_sequence(2);
    ast3 = asteroid_names_sequence(3);
    ast4 = asteroid_names_sequence(4);
    
    %% Position of planet at solutions moments
    [kep_EA, ksun] = uplanet(MJD01, 3); % Earth Departure
    [r1, v1] = sv_from_coe(kep_EA,ksun); %km, km/s
    
    % DSM
    iEA = kep_EA(3);
    rD = sol.rD_mag*[cos(sol.thetaD)*cos(iEA); sin(sol.thetaD)*cos(iEA); sin(iEA)];

    [kep_ast_1_1] = uNEO2(MJDF1,ast1,data);
    [r2,v2] = sv_from_coe(kep_ast_1_1,ksun); % km, km/s
    [kep_ast_1_2] = uNEO2(MJD02,ast1,data);
    [r3,v3] = sv_from_coe(kep_ast_1_2,ksun); % km, km/s % Ast 1 Departure after Coasting
    [kep_ast_2_1] = uNEO2(MJDF2,ast2,data);
    [r4,v4] = sv_from_coe(kep_ast_2_1,ksun); % km, km/s
    [kep_ast_2_2] = uNEO2(MJD03,ast2,data);
    [r5,v5] = sv_from_coe(kep_ast_2_2,ksun); % km, km/s
    [kep_ast_3_1] = uNEO2(MJDF3,ast3,data);
    [r6,v6] = sv_from_coe(kep_ast_3_1,ksun); % km, km/s
    [kep_ast_3_2] = uNEO2(MJD04,ast3,data);
    [r7,v7] = sv_from_coe(kep_ast_3_2,ksun); % km, km/s % Ast 3 Departure after Coasting
    [kep_ast_4_1] = uNEO2(MJDF4,ast4,data);
    [r8,v8] = sv_from_coe(kep_ast_4_1,ksun); % km, km/s % Ast 4 Arrival
    
    % Converting mJ2000 in seconds
t1a_sec = MJD01*60*60*24;
t1b_sec = (MJD01 + sol.TOF1a)*60*60*24;
t2_sec = MJDF1*60*60*24;
t3_sec = MJD02*60*60*24;
t4_sec = MJDF2*60*60*24;
t5_sec = MJD03*60*60*24;
t6_sec = MJDF3*60*60*24;
t7_sec = MJD04*60*60*24;
t8_sec = MJDF4*60*60*24;

% DV calculation with lambert

% Earth -> DSM
[~,~,~,~,VDd,VDa,~,~] = lambertMR(r1,rD,(t1b_sec - t1a_sec),ksun,0,0,0,0);

sol.dv1 = sqrt((VDd(1)- v1(1))^2+(VDd(2)-v1(2))^2+(VDd(3)- v1(3))^2);

if sol.dv1< sqrt(sim.C3_max) % vinf that the launcher can give max 
    sol.dv_extra_launch = 0;
else
    c_launcher = 40; % penalty factor for dv_extra_launch
    sol.dv_extra_launch = c_launcher*(sol.dv1 - sqrt(sim.C3_max))^2; % penalty like, but not discard a priori
end

% DSM -> asteroid 1
[~,~,~,~,V1d,V1a,~,~] = lambertMR(rD,r2,(t2_sec - t1b_sec),ksun,0,0,0,0);
sol.dv2 = sqrt((V1d(1)-VDa(1))^2+(V1d(2)-VDa(2))^2+(V1d(3)-VDa(3))^2);
sol.dv3 = sqrt((V1a(1)- v2(1))^2+(V1a(2)-v2(2))^2+(V1a(3)-v2(3))^2);

% asteroid 1 -> asteroid 2
[~,~,~,~,V2d,V2a,~,~] = lambertMR(r3,r4,(t4_sec - t3_sec),ksun,0,0,0,0);
sol.dv4 = sqrt((V2d(1)-v3(1))^2+(V2d(2)-v3(2))^2+(V2d(3)-v3(3))^2);
sol.dv5 = sqrt((V2a(1)- v4(1))^2+(V2a(2)-v4(2))^2+(V2a(3)-v4(3))^2);

% asteroid 2 -> asteroid 3 
[~,~,~,~,V3d,V3a,~,~] = lambertMR(r5,r6,(t6_sec - t5_sec),ksun,0,0,0,0);
sol.dv6 = sqrt((V3d(1)-v5(1))^2+(V3d(2)-v5(2))^2+(V3d(3)-v5(3))^2);
sol.dv7 = sqrt((V3a(1)- v6(1))^2+(V3a(2)-v6(2))^2+(V3a(3)-v6(3))^2);


% asteroid 3 -> asteroid 4 %%%%%%%%%%
[~,~,~,~,V4d,V4a,~,~] = lambertMR(r7,r8,(t8_sec - t7_sec),ksun,0,0,0,0);
sol.dv8 = sqrt((V4d(1)-v7(1))^2+(V4d(2)-v7(2))^2+(V4d(3)-v7(3))^2);
sol.dv9 = sqrt((V4a(1)- v8(1))^2+(V4a(2)-v8(2))^2+(V4a(3)-v8(3))^2);


sol.dvtot_plot = sol.dv2 + sol.dv3 + sol.dv4+ sol.dv5+ sol.dv6 + sol.dv7 + sol.dv8 + sol.dv9;
    
    %% Plotting
    % PLOT FULL ORBITS AND BEST LAMBERT TRANSFER 
    figure('Name','Mission Orbits and Phases')
    
    % Earth
    plot_earth_orbit(MJD01,3,colors,8);
    hold on
    % Asteroids
    Frac_Orb = 1/6;
    plot_asteorid_orbit(MJDF1,Frac_Orb,ast1,colors,2);
    plot_asteorid_orbit(MJDF2,Frac_Orb,ast2,colors,3);
    plot_asteorid_orbit(MJDF3,Frac_Orb,ast3,colors,4);
    plot_asteorid_orbit(MJDF4,Frac_Orb,ast4,colors,5);
    
    % Mission Arcs
    % First leg: Earth -> DSM
    y012 = [r1; VDd']; %km, km/s; velocity from lambert arc transfer orbit injection
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t12,y12] = ode113(@rates, [t1a_sec t1b_sec], y012,options,'sun');
    plot3( y12(:,1)./AU, y12(:,2)./AU, y12(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','Trajectory');
    
    % 2nd leg: DSM -> AST 1
    y012b = [rD; V1d']; %km, km/s; velocity from lambert arc transfer orbit injection
    [t12b,y12b] = ode113(@rates, [t1b_sec t2_sec], y012b,options,'sun');
    plot3( y12b(:,1)./AU, y12b(:,2)./AU, y12b(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','Trajectory');

    % Coasting 1 on Ast 1
    y023 = [r2; v2]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
    [t23,y23] = ode113(@rates, [t2_sec t3_sec],y023,options,'sun');
    hc1 = plot3( y23(:,1)./AU, y23(:,2)./AU, y23(:,3)./AU,'*','Color',colors(1,:),'Markersize',3);
    hc1.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % 3rd leg: Ast 1 -> Ast 2
    y034 = [r3; V2d']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t34,y34] = ode113(@rates,[t3_sec t4_sec], y034,options,'sun');
    ht2 = plot3( y34(:,1)./AU, y34(:,2)./AU, y34(:,3)./AU,'Color',colors(1,:));
    ht2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    % Coasting 2 on Ast 2
    y045 = [r4; v4]; %km, km/s; 
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t45,y45] = ode113(@rates, [t4_sec t5_sec], y045,options,'sun');
    hc2 = plot3( y45(:,1)./AU, y45(:,2)./AU, y45(:,3)./AU,'*','Color',colors(1,:),'Markersize',3);
    hc2.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % 4th leg: Ast 2 -> Ast 3
    y056 = [r5; V3d']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t56,y56] = ode113(@rates,[t5_sec t6_sec], y056,options,'sun');
    ht3 = plot3( y56(:,1)./AU, y56(:,2)./AU, y56(:,3)./AU,'Color',colors(1,:));
    ht3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    % Coasting 3 on Ast 3
    y067 = [r6; v6]; %km, km/s; after the application of 2nd dv of its transfer leg, rendezvous with asteroid, same velocity as the asteroid at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t67,y67] = ode113(@rates,[t6_sec t7_sec], y067,options,'sun');
    hc3 = plot3( y67(:,1)./AU, y67(:,2)./AU, y67(:,3)./AU,'*','Color',colors(1,:),'Markersize',3);
    hc3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    % 5th leg: Ast 3 -> Ast 4
    y078 = [r7; V4d']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t78,y78] = ode113(@rates, [t7_sec t8_sec], y078,options,'sun');
    ht4 = plot3( y78(:,1)./AU, y78(:,2)./AU, y78(:,3)./AU,'Color',colors(1,:));
    ht4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    

    %% Naming
    % Sun Yellow Asterisk
    plot3(0,0,0,'*','Color',colors(4,:),'DisplayName','Sun');
    
    legend('show','Location','southeastoutside')
    
    hp1 = plot3(r1(1)./AU,r1(2)./AU,r1(3)./AU,'o','Color',colors(8,:),'MarkerSize',6,...
    'DisplayName','Earth Departure');

    hp1b = plot3(rD(1)./AU,rD(2)./AU,rD(3)./AU,'o','Color',colors(12,:),'MarkerSize',6,...
    'DisplayName','DSM');

    name_hp2 = sprintf('Arrival %s',ast1);
    hp2 = plot3(r2(1)./AU,r2(2)./AU,r2(3)./AU,'^','Color',colors(2,:),'MarkerSize',6,...
    'DisplayName',name_hp2);
%     hp2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp3 = sprintf('Departure %s',ast1);
    hp3 = plot3(r3(1)./AU,r3(2)./AU,r3(3)./AU,'o','Color',colors(2,:),'MarkerSize',6,...
    'DisplayName',name_hp3);
%     hp3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp4 = sprintf('Arrival %s',ast2);
    hp4 = plot3(r4(1)./AU,r4(2)./AU,r4(3)./AU,'^','Color',colors(3,:),'MarkerSize',6,...
    'DisplayName',name_hp4);
%     hp4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp5 = sprintf('Departure %s',ast2);
    hp5 = plot3(r5(1)./AU,r5(2)./AU,r5(3)./AU,'o','Color',colors(3,:),'MarkerSize',6,...
    'DisplayName',name_hp5);
%     hp5.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp6 = sprintf('Arrival %s',ast3);
    hp6 = plot3(r6(1)./AU,r6(2)./AU,r6(3)./AU,'^','Color',colors(4,:),'MarkerSize',6,...
    'DisplayName',name_hp6);
%     hp6.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp7 = sprintf('Departure %s',ast3);
    hp7 = plot3(r7(1)./AU,r7(2)./AU,r7(3)./AU,'o','Color',colors(4,:),'MarkerSize',6,...
    'DisplayName',name_hp7);
%     hp7.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp8 = sprintf('Arrival %s',ast4);
    hp8 = plot3(r8(1)./AU,r8(2)./AU,r8(3)./AU,'^','Color',colors(5,:),'MarkerSize',6,...
    'DisplayName',name_hp8);
%     hp8.Annotation.LegendInformation.IconDisplayStyle = 'off';

    axis equal; grid on
    title(sim.case_name)
    xlabel('x [AU]')
    ylabel('y [AU]')
    zlabel('z [AU]')

end