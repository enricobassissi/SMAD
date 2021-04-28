function [sol] = plot_mission_4neo(sol,asteroid_names_sequence,data,sim,colors)
    
    AU = astroConstants(2);

    % initialise stuff from solution
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
    
    % Position of planet at solutions moments
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
    
    % Converting mJD2000 in seconds
    t1_sec = MJD01*60*60*24;
    t2_sec = MJDF1*60*60*24;
    t3_sec = MJD02*60*60*24;
    t4_sec = MJDF2*60*60*24;
    t5_sec = MJD03*60*60*24;
    t6_sec = MJDF3*60*60*24;
    t7_sec = MJD04*60*60*24;
    t8_sec = MJDF4*60*60*24;
    
%     % Earth injection orbits
%     v_inf = sol.v_inf_magn.*[cos(sol.v_inf_alpha)*cos(sol.v_inf_beta); 
%                     sin(sol.v_inf_alpha)*cos(sol.v_inf_beta); 
%                     sin(sol.v_inf_beta)];
%     v_inj = v_EA+v_inf'; % v_EA is row, v_inf would be column as built
    
    % Lamberts and deltaVs
    ToF12_sec = t2_sec - t1_sec;
    [~,~,~,~,VI12,VF12,~,~] = lambertMR(r_EA,r_ast1_1,ToF12_sec,ksun,0,0,0,0);
    dv1 = sqrt((VI12(1)-v_EA(1))^2+(VI12(2)-v_EA(2))^2+(VI12(3)-v_EA(3))^2);
    if dv1 < sqrt(sim.C3_max)
        sol.dV_single.dV_extra_launch = 0;
        sol.Vinf_launcher = dv1;
    else
        sol.dV_single.dV_extra_launch = dv1 - sqrt(sim.C3_max);
        sol.Vinf_launcher = sqrt(sim.C3_max);
    end
    sol.dV_single.dV2 = sqrt((VF12(1)-v_ast1_1(1))^2+(VF12(2)-v_ast1_1(2))^2+(VF12(3)-v_ast1_1(3))^2);
    sol.dV_tot_leg1 = sol.dV_single.dV2 + sol.dV_single.dV_extra_launch;
    
    ToF34_sec = t4_sec - t3_sec;
    [~,~,~,~,VI34,VF34,~,~] = lambertMR(r_ast1_2,r_ast2_1,ToF34_sec,ksun,0,0,0,0);
    sol.dV_single.dV3 = sqrt((VI34(1)-v_ast1_2(1))^2+(VI34(2)-v_ast1_2(2))^2+(VI34(3)-v_ast1_2(3))^2);
    sol.dV_single.dV4 = sqrt((VF34(1)-v_ast2_1(1))^2+(VF34(2)-v_ast2_1(2))^2+(VF34(3)-v_ast2_1(3))^2);
    sol.dV_tot_leg2 = sol.dV_single.dV3 + sol.dV_single.dV4;
    
    ToF56_sec = t6_sec - t5_sec;
    [~,~,~,~,VI56,VF56,~,~] = lambertMR(r_ast2_2,r_ast3_1,ToF56_sec,ksun,0,0,0,0);
    sol.dV_single.dV5 = sqrt((VI56(1)-v_ast2_2(1))^2+(VI56(2)-v_ast2_2(2))^2+(VI56(3)-v_ast2_2(3))^2);
    sol.dV_single.dV6 = sqrt((VF56(1)-v_ast3_1(1))^2+(VF56(2)-v_ast3_1(2))^2+(VF56(3)-v_ast3_1(3))^2);
    sol.dV_tot_leg3 = sol.dV_single.dV5 + sol.dV_single.dV6;
    
    ToF78_sec = t8_sec - t7_sec;
    [~,~,~,~,VI78,VF78,~,~] = lambertMR(r_ast3_2,r_ast4_1,ToF78_sec,ksun,0,0,0,0);
    sol.dV_single.dV7 = sqrt((VI78(1)-v_ast3_2(1))^2+(VI78(2)-v_ast3_2(2))^2+(VI78(3)-v_ast3_2(3))^2);
    sol.dV_single.dV8 = sqrt((VF78(1)-v_ast4_1(1))^2+(VF78(2)-v_ast4_1(2))^2+(VF78(3)-v_ast4_1(3))^2);
    sol.dV_tot_leg4 = sol.dV_single.dV7 + sol.dV_single.dV8;

    % PLOT FULL ORBITS AND BEST LAMBERT TRANSFER 
    figure('Name','Mission Orbits and Phases')
    % Earth
    plot_earth_orbit(MJD01,colors,8);
    hold on
    % Asteroids
    years = 6;
    plot_asteorid_orbit(MJDF1,years,ast1,colors,2);
    plot_asteorid_orbit(MJDF2,years,ast2,colors,3);
    plot_asteorid_orbit(MJDF3,years,ast3,colors,4);
    plot_asteorid_orbit(MJDF4,years,ast4,colors,5);
    
    % Mission Arcs
    % First leg: Earth -> Ast 1
    t012 = t1_sec;
    tf12 = t2_sec;
    y012 = [r_EA; VI12']; %km, km/s; velocity from lambert arc transfer orbit injection
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y12] = ode113(@rates, [t012 tf12], y012,options,'sun');
    plot3( y12(:,1)./AU, y12(:,2)./AU, y12(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','First Leg');

    % Coasting 1 on Ast 1
    t023 = t2_sec;
    tf23 = t3_sec;
    y023 = [r_ast1_1; v_ast1_1]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y23] = ode113(@rates, [t023 tf23],y023,options,'sun');
    plot3( y23(:,1)./AU, y23(:,2)./AU, y23(:,3)./AU,'*','Color',colors(1,:),'Markersize',3,...
        'DisplayName','Coasting 1');

    % Second leg: Ast 1 -> Ast 2
    t034 = t3_sec;
    tf34 = t4_sec;
    y034 = [r_ast1_2; VI34']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y34] = ode113(@rates, [t034 tf34], y034,options,'sun');
    plot3( y34(:,1)./AU, y34(:,2)./AU, y34(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','Leg 2');
    
    % Coasting 2 on Ast 2
    t045 = t4_sec;
    tf45 = t5_sec;
    y045 = [r_ast2_1; v_ast2_1]; %km, km/s; after the application of 2nd dv of its transfer leg, rendezvous with asteroid, same velocity as the asteroid at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y45] = ode113(@rates, [t045 tf45], y045,options,'sun');
    plot3( y45(:,1)./AU, y45(:,2)./AU, y45(:,3)./AU,'*','Color',colors(1,:),'Markersize',3,...
        'DisplayName','Coasting 2');

    % Third leg: Ast 2 -> Ast 3
    t056 = t5_sec;
    tf56 = t6_sec;
    y056 = [r_ast2_2; VI56']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y56] = ode113(@rates, [t056 tf56], y056,options,'sun');
    plot3( y56(:,1)./AU, y56(:,2)./AU, y56(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','Leg 3');
    
    % Coasting 3 on Ast 3
    t067 = t6_sec;
    tf67 = t7_sec;
    y067 = [r_ast3_1; v_ast3_1]; %km, km/s; after the application of 2nd dv of its transfer leg, rendezvous with asteroid, same velocity as the asteroid at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y67] = ode113(@rates, [t067 tf67], y067,options,'sun');
    plot3( y67(:,1)./AU, y67(:,2)./AU, y67(:,3)./AU,'*','Color',colors(1,:),'Markersize',3,...
        'DisplayName','Coasting 3');

    % Fourth leg: Ast 3 -> Ast 4
    t078 = t7_sec;
    tf78 = t8_sec;
    y078 = [r_ast3_2; VI78']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y78] = ode113(@rates, [t078 tf78], y078,options,'sun');
    plot3( y78(:,1)./AU, y78(:,2)./AU, y78(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','Leg 4');
    
    % Sun Yellow Asterisk
    plot3(0,0,0,'*','Color',colors(4,:),'DisplayName','Sun');
    
    legend('show','Location','southeastoutside')
    
    hp1 = plot3(r_EA(1)./AU,r_EA(2)./AU,r_EA(3)./AU,'o','Color',colors(8,:),'MarkerSize',4);
    hp1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp2 = plot3(r_ast1_1(1)./AU,r_ast1_1(2)./AU,r_ast1_1(3)./AU,'^','Color',colors(2,:),'MarkerSize',4);
    hp2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp3 = plot3(r_ast1_2(1)./AU,r_ast1_2(2)./AU,r_ast1_2(3)./AU,'o','Color',colors(2,:),'MarkerSize',4);
    hp3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp4 = plot3(r_ast2_1(1)./AU,r_ast2_1(2)./AU,r_ast2_1(3)./AU,'^','Color',colors(3,:),'MarkerSize',4);
    hp4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp5 = plot3(r_ast2_2(1)./AU,r_ast2_2(2)./AU,r_ast2_2(3)./AU,'o','Color',colors(3,:),'MarkerSize',4);
    hp5.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp6 = plot3(r_ast3_1(1)./AU,r_ast3_1(2)./AU,r_ast3_1(3)./AU,'^','Color',colors(4,:),'MarkerSize',4);
    hp6.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp7 = plot3(r_ast3_2(1)./AU,r_ast3_2(2)./AU,r_ast3_2(3)./AU,'o','Color',colors(4,:),'MarkerSize',4);
    hp7.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp8 = plot3(r_ast4_1(1)./AU,r_ast4_1(2)./AU,r_ast4_1(3)./AU,'^','Color',colors(5,:),'MarkerSize',4);
    hp8.Annotation.LegendInformation.IconDisplayStyle = 'off';

    axis equal; grid on
    xlabel('AU')
    ylabel('AU')
    zlabel('AU')

end