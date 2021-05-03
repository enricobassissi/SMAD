function [sol] = plot_mission_4neo_flyby(sol,asteroid_names_sequence,data,sim,colors)
    
    AU = astroConstants(2);
    
    %% Initialise stuff from solution
    MJD01 = sol.MJD0;
    MJDP1 = MJD01 + sol.TOF1;
    MJDP2 = MJDP1 + sol.TOF2;
    MJDP3 = MJDP2 + sol.TOF3;
    MJDP4 = MJDP3 + sol.TOF4;
    
    ast1 = asteroid_names_sequence(1);
    ast2 = asteroid_names_sequence(2);
    ast3 = asteroid_names_sequence(3);
    ast4 = asteroid_names_sequence(4);
    
    %% Position of planet at solutions moments
    [kep_EA, ksun] = uplanet(MJD01, 3); % Earth Departure
    [rEA, vEA] = sv_from_coe(kep_EA,ksun); %km, km/s
    [kep_ast_1] = uNEO2(MJDP1,ast1,data);
    [r_ast1,v_ast1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
    [kep_ast_2] = uNEO2(MJDP2,ast2,data);
    [r_ast2,v_ast2] = sv_from_coe(kep_ast_2,ksun); % km, km/s
    [kep_ast_3] = uNEO2(MJDP3,ast3,data);
    [r_ast3,v_ast3] = sv_from_coe(kep_ast_3,ksun); % km, km/s
    [kep_ast_4] = uNEO2(MJDP4,ast4,data);
    [r_ast4,v_ast4] = sv_from_coe(kep_ast_4,ksun); % km, km/s % Ast 4 Arrival
    
    % Converting mJD2000 in seconds
    t1_sec = MJD01*60*60*24;
    t2_sec = MJDP1*60*60*24;
    t3_sec = MJDP2*60*60*24;
    t4_sec = MJDP3*60*60*24;
    t5_sec = MJDP4*60*60*24;
    
    %% Lamberts and deltaVs
    ToF_EAast1_sec = t2_sec - t1_sec;
    [~,~,~,~,VI_EAast1,VF_EAast1,~,~] = lambertMR(rEA,r_ast1,ToF_EAast1_sec,ksun,0,0,0,0);
    dv1 = sqrt((VI_EAast1(1)-vEA(1))^2+(VI_EAast1(2)-vEA(2))^2+(VI_EAast1(3)-vEA(3))^2);
    if dv1 < sqrt(sim.C3_max)
        sol.dV_extra_launch = 0;
        sol.Vinf_launcher = dv1;
    else
        sol.dV_extra_launch = dv1 - sqrt(sim.C3_max);
        sol.Vinf_launcher = sqrt(sim.C3_max);
    end
%     sol.dV_tot_leg1 = sol.dV_single.dV2 + sol.dV_single.dV_extra_launch;
    % relative velocity arrival at the asteroid and asteroid itself, for the deployment of the "lander"
    sol.Vrel_passage_ast1 = sqrt((VF_EAast1(1)-v_ast1(1))^2+(VF_EAast1(2)-v_ast1(2))^2+(VF_EAast1(3)-v_ast1(3))^2);
    
    ToF_ast12_sec = t3_sec - t2_sec;
    [~,~,~,~,VI_ast12,VF_ast12,~,~] = lambertMR(r_ast1,r_ast2,ToF_ast12_sec,ksun,0,0,0,0);
%     sol.dV_single.dV3 = sqrt((VI34(1)-v_ast1_2(1))^2+(VI34(2)-v_ast1_2(2))^2+(VI34(3)-v_ast1_2(3))^2);
%     sol.dV_single.dV4 = sqrt((VF34(1)-v_ast2_1(1))^2+(VF34(2)-v_ast2_1(2))^2+(VF34(3)-v_ast2_1(3))^2);
%     sol.dV_tot_leg2 = sol.dV_single.dV3 + sol.dV_single.dV4;
    
    % dV of flyby passage on asteroid 1
    sol.dVast1 = sqrt((VI_ast12(1)-VF_EAast1(1))^2+(VI_ast12(2)-VF_EAast1(2))^2+(VI_ast12(3)-VF_EAast1(3))^2);
    % relative velocity arrival at the asteroid and asteroid itself, for the deployment of the "lander"
    sol.Vrel_passage_ast2 = sqrt((VF_ast12(1)-v_ast2(1))^2+(VF_ast12(2)-v_ast2(2))^2+(VF_ast12(3)-v_ast2(3))^2);
    
    ToF_ast23_sec = t4_sec - t3_sec;
    [~,~,~,~,VI_ast23,VF_ast23,~,~] = lambertMR(r_ast2,r_ast3,ToF_ast23_sec,ksun,0,0,0,0);
%     sol.dV_single.dV5 = sqrt((VI56(1)-v_ast2_2(1))^2+(VI56(2)-v_ast2_2(2))^2+(VI56(3)-v_ast2_2(3))^2);
%     sol.dV_tot_leg3 = sol.dV_single.dV5 + sol.dV_single.dV6;

    % dV of flyby passage on asteroid 2
    sol.dVast2 = sqrt((VI_ast23(1)-VF_ast12(1))^2+(VI_ast23(2)-VF_ast12(2))^2+(VI_ast23(3)-VF_ast12(3))^2);
    % relative velocity arrival at the asteroid and asteroid itself, for the deployment of the "lander"
    sol.Vrel_passage_ast3 = sqrt((VF_ast23(1)-v_ast3(1))^2+(VF_ast23(2)-v_ast3(2))^2+(VF_ast23(3)-v_ast3(3))^2);
    
    ToF_ast34_sec = t5_sec - t4_sec;
    [~,~,~,~,VI_ast34,VF_ast34,~,~] = lambertMR(r_ast3,r_ast4,ToF_ast34_sec,ksun,0,0,0,0);
%     sol.dV_single.dV7 = sqrt((VI78(1)-v_ast3_2(1))^2+(VI78(2)-v_ast3_2(2))^2+(VI78(3)-v_ast3_2(3))^2);
%     sol.dV_single.dV8 = sqrt((VF78(1)-v_ast4_1(1))^2+(VF78(2)-v_ast4_1(2))^2+(VF78(3)-v_ast4_1(3))^2);
%     sol.dV_tot_leg4 = sol.dV_single.dV7 + sol.dV_single.dV8;

    % dV of flyby passage on asteroid 3
    sol.dVast3 = sqrt((VI_ast34(1)-VF_ast23(1))^2+(VI_ast34(2)-VF_ast23(2))^2+(VI_ast34(3)-VF_ast23(3))^2);
    % relative velocity arrival at the asteroid and asteroid itself, for the deployment of the "lander"
    sol.Vrel_passage_ast4 = sqrt((VF_ast34(1)-v_ast4(1))^2+(VF_ast34(2)-v_ast4(2))^2+(VF_ast34(3)-v_ast4(3))^2);
    
    % actual total dV, different from Fval because now it enters also
    % penalty in the optimisation function
    sol.dV_tot = sol.dV_extra_launch + sol.dVast1 + sol.dVast2 + sol.dVast3;
    
    %% Plotting
    % PLOT FULL ORBITS AND BEST LAMBERT TRANSFER 
    figure('Name','Mission Orbits and Phases')
    % Earth
    plot_earth_orbit(MJD01,colors,8);
    hold on
    % Asteroids
    years = 5;
%     plot_asteorid_orbit(MJDP1,years,ast1,colors,2);
%     plot_asteorid_orbit(MJDP2,years,ast2,colors,3);
%     plot_asteorid_orbit(MJDP3,years,ast3,colors,4);
%     plot_asteorid_orbit(MJDP4,years,ast4,colors,5);
    
    % Mission Arcs
    % First leg: Earth -> Ast 1
    t012 = t1_sec;
    tf12 = t2_sec;
    y012 = [rEA; VI_EAast1']; %km, km/s; velocity from lambert arc transfer orbit injection
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t12,y12] = ode113(@rates, [t012 tf12], y012,options,'sun');
    plot3( y12(:,1)./AU, y12(:,2)./AU, y12(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','Trajectory');

    % Second leg: Ast 1 -> Ast 2
    t034 = t2_sec;
    tf34 = t3_sec;
    y034 = [r_ast1; VI_ast12']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t34,y34] = ode113(@rates, [t034 tf34], y034,options,'sun');
    ht2 = plot3( y34(:,1)./AU, y34(:,2)./AU, y34(:,3)./AU,'Color',colors(1,:));
    ht2.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % Third leg: Ast 2 -> Ast 3
    t056 = t3_sec;
    tf56 = t4_sec;
    y056 = [r_ast2; VI_ast23']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t56,y56] = ode113(@rates, [t056 tf56], y056,options,'sun');
    ht3 = plot3( y56(:,1)./AU, y56(:,2)./AU, y56(:,3)./AU,'Color',colors(1,:));
    ht3.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % Fourth leg: Ast 3 -> Ast 4
    t078 = t4_sec;
    tf78 = t5_sec;
    y078 = [r_ast3; VI_ast34']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t78,y78] = ode113(@rates, [t078 tf78], y078,options,'sun');
    ht4 =plot3( y78(:,1)./AU, y78(:,2)./AU, y78(:,3)./AU,'Color',colors(1,:));
    ht4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    %% Extract the Sun-SpaceCraft Distance for all the trajectory [km]
    sol.SunSpacecraftDistanceNorm = vecnorm([y12;y34;y56;y78],2,2); % 2,2 means norm 2 and by row
    sol.SpacecraftTrajectory = [y12;y34;y56;y78];
    sol.SCtime = [t12;t34;t56;t78];
    [sol.angles.SAA,sol.angles.EVA,sol.angles.SCA,sol.angles.SolarConjunction] = aspect_angles(sol);
    
    %% Naming
    % Sun Yellow Asterisk
    plot3(0,0,0,'*','Color',colors(4,:),'DisplayName','Sun');
    
    legend('show','Location','southeastoutside')
    
    hp1 = plot3(rEA(1)./AU,rEA(2)./AU,rEA(3)./AU,'o','Color',colors(8,:),'MarkerSize',6,...
    'DisplayName','Earth Departure');
%     hp1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp2 = plot3(r_ast1(1)./AU,r_ast1(2)./AU,r_ast1(3)./AU,'^','Color',colors(2,:),'MarkerSize',6,...
        'DisplayName',ast1);
%     hp2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp4 = plot3(r_ast2(1)./AU,r_ast2(2)./AU,r_ast2(3)./AU,'^','Color',colors(3,:),'MarkerSize',6,...
        'DisplayName',ast2);
%     hp4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp6 = plot3(r_ast3(1)./AU,r_ast3(2)./AU,r_ast3(3)./AU,'^','Color',colors(4,:),'MarkerSize',6,...
        'DisplayName',ast3);
%     hp6.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp8 = plot3(r_ast4(1)./AU,r_ast4(2)./AU,r_ast4(3)./AU,'^','Color',colors(5,:),'MarkerSize',6,...
        'DisplayName',ast4);
%     hp8.Annotation.LegendInformation.IconDisplayStyle = 'off';

    axis equal; grid on
    title(sim.case_name)
    xlabel('AU')
    ylabel('AU')
    zlabel('AU')
    
end