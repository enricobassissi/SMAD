function [sol] = plot_mission_4neo_rendezvous(sol,asteroid_names_sequence,data,sim,colors)
    
    AU = astroConstants(2);

    %% initialise stuff from solution
    MJD01 = sol.MJD0;
    MJDF1 = MJD01 + sol.TOF1;
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
    [r_EA, v_EA] = sv_from_coe(kep_EA,ksun); %km, km/s
    [kep_ast_1_1] = uNEO2(MJDF1,ast1,data);
    [r_ast1_1,v_ast1_1] = sv_from_coe(kep_ast_1_1,ksun); % km, km/s
    [kep_ast_1_2] = uNEO2(MJD02,ast1,data);
    [r_ast1_2,v_ast1_2] = sv_from_coe(kep_ast_1_2,ksun); % km, km/s % Ast 1 Departure after Coasting
    [kep_ast_2_1] = uNEO2(MJDF2,ast2,data);
    [r_ast2_1,v_ast2_1] = sv_from_coe(kep_ast_2_1,ksun); % km, km/s
    [kep_ast_2_2] = uNEO2(MJD03,ast2,data);
    [r_ast2_2,v_ast2_2] = sv_from_coe(kep_ast_2_2,ksun); % km, km/s
    [kep_ast_3_1] = uNEO2(MJDF3,ast3,data);
    [r_ast3_1,v_ast3_1] = sv_from_coe(kep_ast_3_1,ksun); % km, km/s
    [kep_ast_3_2] = uNEO2(MJD04,ast3,data);
    [r_ast3_2,v_ast3_2] = sv_from_coe(kep_ast_3_2,ksun); % km, km/s % Ast 3 Departure after Coasting
    [kep_ast_4_1] = uNEO2(MJDF4,ast4,data);
    [r_ast4_1,v_ast4_1] = sv_from_coe(kep_ast_4_1,ksun); % km, km/s % Ast 4 Arrival
    
    % Converting time points from days to seconds
    departure_sec = sol.MJD0*60*60*24;
    ToF12_sec = sol.TOF1*60*60*24;
    buffer_time1_sec = sol.buffer_time1*60*60*24;
    ToF34_sec = sol.TOF2*60*60*24;
    buffer_time2_sec = sol.buffer_time2*60*60*24;
    ToF56_sec = sol.TOF3*60*60*24;
    buffer_time3_sec = sol.buffer_time3*60*60*24;
    ToF78_sec = sol.TOF4*60*60*24;
    
    %% Lamberts and deltaVs
    [~,~,~,~,VI12,VF12,~,~] = lambertMR(r_EA,r_ast1_1,ToF12_sec,ksun,0,0,0,0);
    dv1 = sqrt((VI12(1)-v_EA(1))^2+(VI12(2)-v_EA(2))^2+(VI12(3)-v_EA(3))^2);
    if dv1 < sqrt(sim.C3_max)
        sol.dV_extra_launch = 0;
        sol.Vinf_launcher = dv1;
    else
        sol.dV_extra_launch = dv1 - sqrt(sim.C3_max);
        sol.Vinf_launcher = sqrt(sim.C3_max);
    end
    sol.dV2 = sqrt((VF12(1)-v_ast1_1(1))^2+(VF12(2)-v_ast1_1(2))^2+(VF12(3)-v_ast1_1(3))^2);
    sol.dV_tot_leg1 = sol.dV2 + sol.dV_extra_launch;
    
    [~,~,~,~,VI34,VF34,~,~] = lambertMR(r_ast1_2,r_ast2_1,ToF34_sec,ksun,0,0,0,0);
    sol.dV3 = sqrt((VI34(1)-v_ast1_2(1))^2+(VI34(2)-v_ast1_2(2))^2+(VI34(3)-v_ast1_2(3))^2);
    sol.dV4 = sqrt((VF34(1)-v_ast2_1(1))^2+(VF34(2)-v_ast2_1(2))^2+(VF34(3)-v_ast2_1(3))^2);
    sol.dV_tot_leg2 = sol.dV3 + sol.dV4;
    
    [~,~,~,~,VI56,VF56,~,~] = lambertMR(r_ast2_2,r_ast3_1,ToF56_sec,ksun,0,0,0,0);
    sol.dV5 = sqrt((VI56(1)-v_ast2_2(1))^2+(VI56(2)-v_ast2_2(2))^2+(VI56(3)-v_ast2_2(3))^2);
    sol.dV6 = sqrt((VF56(1)-v_ast3_1(1))^2+(VF56(2)-v_ast3_1(2))^2+(VF56(3)-v_ast3_1(3))^2);
    sol.dV_tot_leg3 = sol.dV5 + sol.dV6;
    
    [~,~,~,~,VI78,VF78,~,~] = lambertMR(r_ast3_2,r_ast4_1,ToF78_sec,ksun,0,0,0,0);
    sol.dV7 = sqrt((VI78(1)-v_ast3_2(1))^2+(VI78(2)-v_ast3_2(2))^2+(VI78(3)-v_ast3_2(3))^2);
    sol.dV8 = sqrt((VF78(1)-v_ast4_1(1))^2+(VF78(2)-v_ast4_1(2))^2+(VF78(3)-v_ast4_1(3))^2);
    sol.dV_tot_leg4 = sol.dV7 + sol.dV8;

    sol.dV_tot = sol.dV_tot_leg1+sol.dV_tot_leg2+sol.dV_tot_leg3+sol.dV_tot_leg4;
    
    %% Plotting
    % PLOT FULL ORBITS AND BEST LAMBERT TRANSFER 
    figure('Name','Mission Orbits and Phases')
    
    % Earth
    plot_earth_orbit(MJD01,colors,8);
    hold on
    % Asteroids
    years = 6;
%     plot_asteorid_orbit(MJDF1,years,ast1,colors,2);
%     plot_asteorid_orbit(MJDF2,years,ast2,colors,3);
%     plot_asteorid_orbit(MJDF3,years,ast3,colors,4);
%     plot_asteorid_orbit(MJDF4,years,ast4,colors,5);
    
    % Mission Arcs
    % First leg: Earth -> Ast 1
    departure_EA_sec = departure_sec;
    arrival_ast1_sec = departure_EA_sec+ToF12_sec;
    time_int12 = [departure_EA_sec, arrival_ast1_sec];
    y012 = [r_EA; VI12']; %km, km/s; velocity from lambert arc transfer orbit injection
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t12,y12] = ode113(@rates, time_int12, y012,options,'sun');
    plot3( y12(:,1)./AU, y12(:,2)./AU, y12(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','Trajectory');

    % Coasting 1 on Ast 1
    start_ast1_sec = arrival_ast1_sec;
    end_ast1_sec = start_ast1_sec+buffer_time1_sec;
    time_int23 = [start_ast1_sec, end_ast1_sec];
    y023 = [r_ast1_1; v_ast1_1]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t23,y23] = ode113(@rates, time_int23,y023,options,'sun');
    hc1 = plot3( y23(:,1)./AU, y23(:,2)./AU, y23(:,3)./AU,'*','Color',colors(1,:),'Markersize',3);
    hc1.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % Second leg: Ast 1 -> Ast 2
    departure_ast1_sec = end_ast1_sec;
    arrival_ast2_sec = departure_ast1_sec+ToF34_sec;
    time_int34 = [departure_ast1_sec, arrival_ast2_sec];
    y034 = [r_ast1_2; VI34']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t34,y34] = ode113(@rates,time_int34, y034,options,'sun');
    ht2 = plot3( y34(:,1)./AU, y34(:,2)./AU, y34(:,3)./AU,'Color',colors(1,:));
    ht2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    % Coasting 2 on Ast 2
    start_ast2_sec = arrival_ast2_sec;
    end_ast2_sec = start_ast2_sec+buffer_time2_sec;
    time_int45 = [start_ast2_sec, end_ast2_sec];
    y045 = [r_ast2_1; v_ast2_1]; %km, km/s; after the application of 2nd dv of its transfer leg, rendezvous with asteroid, same velocity as the asteroid at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t45,y45] = ode113(@rates, time_int45, y045,options,'sun');
    hc2 = plot3( y45(:,1)./AU, y45(:,2)./AU, y45(:,3)./AU,'*','Color',colors(1,:),'Markersize',3);
    hc2.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % Third leg: Ast 2 -> Ast 3
    departure_ast2_sec = end_ast2_sec;
    arrival_ast3_sec = departure_ast2_sec+ToF56_sec;
    time_int56 = [departure_ast2_sec, arrival_ast3_sec];
    y056 = [r_ast2_2; VI56']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t56,y56] = ode113(@rates,time_int56, y056,options,'sun');
    ht3 = plot3( y56(:,1)./AU, y56(:,2)./AU, y56(:,3)./AU,'Color',colors(1,:));
    ht3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    % Coasting 3 on Ast 3
    start_ast3_sec = arrival_ast3_sec;
    end_ast3_sec = start_ast3_sec+buffer_time3_sec;
    time_int67 = [start_ast3_sec, end_ast3_sec];
    y067 = [r_ast3_1; v_ast3_1]; %km, km/s; after the application of 2nd dv of its transfer leg, rendezvous with asteroid, same velocity as the asteroid at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t67,y67] = ode113(@rates,time_int67, y067,options,'sun');
    hc3 = plot3( y67(:,1)./AU, y67(:,2)./AU, y67(:,3)./AU,'*','Color',colors(1,:),'Markersize',3);
    hc3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    % Fourth leg: Ast 3 -> Ast 4
    departure_ast3_sec = end_ast3_sec;
    arrival_ast4_sec = departure_ast3_sec+ToF78_sec;
    time_int78 = [departure_ast3_sec, arrival_ast4_sec];
    y078 = [r_ast3_2; VI78']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t78,y78] = ode113(@rates, time_int78, y078,options,'sun');
    ht4 = plot3( y78(:,1)./AU, y78(:,2)./AU, y78(:,3)./AU,'Color',colors(1,:));
    ht4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    %% Extract the Sun-SpaceCraft Distance for all the trajectory [km]
    sol.SunSpacecraftDistanceNorm = vecnorm([y12;y23;y34;y45;y56;y67;y78],2,2); % 2,2 means norm 2 and by row
    sol.SpacecraftTrajectory = [y12;y23;y34;y45;y56;y67;y78];
    sol.SCtime = [t12;t23;t34;t45;t56;t67;t78];
    [sol.angles.SAA,sol.angles.EVA,sol.angles.SCA,sol.angles.SolarConjunction] = aspect_angles(sol);
    
    %% Naming
    % Sun Yellow Asterisk
    plot3(0,0,0,'*','Color',colors(4,:),'DisplayName','Sun');
    
    legend('show','Location','southeastoutside')
    
    hp1 = plot3(r_EA(1)./AU,r_EA(2)./AU,r_EA(3)./AU,'o','Color',colors(8,:),'MarkerSize',6,...
    'DisplayName','Earth Departure');
%     hp1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp2 = sprintf('Arrival %s',ast1);
    hp2 = plot3(r_ast1_1(1)./AU,r_ast1_1(2)./AU,r_ast1_1(3)./AU,'^','Color',colors(2,:),'MarkerSize',6,...
    'DisplayName',name_hp2);
%     hp2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp3 = sprintf('Departure %s',ast1);
    hp3 = plot3(r_ast1_2(1)./AU,r_ast1_2(2)./AU,r_ast1_2(3)./AU,'o','Color',colors(2,:),'MarkerSize',6,...
    'DisplayName',name_hp3);
%     hp3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp4 = sprintf('Arrival %s',ast2);
    hp4 = plot3(r_ast2_1(1)./AU,r_ast2_1(2)./AU,r_ast2_1(3)./AU,'^','Color',colors(3,:),'MarkerSize',6,...
    'DisplayName',name_hp4);
%     hp4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp5 = sprintf('Departure %s',ast2);
    hp5 = plot3(r_ast2_2(1)./AU,r_ast2_2(2)./AU,r_ast2_2(3)./AU,'o','Color',colors(3,:),'MarkerSize',6,...
    'DisplayName',name_hp5);
%     hp5.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp6 = sprintf('Arrival %s',ast3);
    hp6 = plot3(r_ast3_1(1)./AU,r_ast3_1(2)./AU,r_ast3_1(3)./AU,'^','Color',colors(4,:),'MarkerSize',6,...
    'DisplayName',name_hp6);
%     hp6.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp7 = sprintf('Departure %s',ast3);
    hp7 = plot3(r_ast3_2(1)./AU,r_ast3_2(2)./AU,r_ast3_2(3)./AU,'o','Color',colors(4,:),'MarkerSize',6,...
    'DisplayName',name_hp7);
%     hp7.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp8 = sprintf('Arrival %s',ast4);
    hp8 = plot3(r_ast4_1(1)./AU,r_ast4_1(2)./AU,r_ast4_1(3)./AU,'^','Color',colors(5,:),'MarkerSize',6,...
    'DisplayName',name_hp8);
%     hp8.Annotation.LegendInformation.IconDisplayStyle = 'off';

    axis equal; grid on
    title(sim.case_name)
    xlabel('AU')
    ylabel('AU')
    zlabel('AU')

end