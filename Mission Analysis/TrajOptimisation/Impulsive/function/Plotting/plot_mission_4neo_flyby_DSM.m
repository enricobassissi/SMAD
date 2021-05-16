function [sol] = plot_mission_4neo_flyby_DSM(sol,asteroid_names_sequence,data,sim,colors)
    
    AU = astroConstants(2);
    
    %% Initialise stuff from solution
    MJD01 = sol.MJD0;
    MJDP1a = MJD01 + sol.TOF1a;
    MJDP1b = MJDP1a + sol.TOF1b;
    MJDP2 = MJDP1b + sol.TOF2;
    MJDP3 = MJDP2 + sol.TOF3;
    MJDP4 = MJDP3 + sol.TOF4;
    
    ast1 = asteroid_names_sequence(1);
    ast2 = asteroid_names_sequence(2);
    ast3 = asteroid_names_sequence(3);
    ast4 = asteroid_names_sequence(4);
    
    %% Position of planet at solutions moments
    [kep_EA, ksun] = uplanet(MJD01, 3); % Earth Departure
    [rEA, vEA] = sv_from_coe(kep_EA,ksun); %km, km/s
    
    % DSM
    iEA = kep_EA(3);
    rD = sol.rD_mag*[cos(sol.thetaD)*cos(iEA); sin(sol.thetaD)*cos(iEA); sin(iEA)];
    
    [kep_ast_1] = uNEO2(MJDP1b,ast1,data);
    [r1,v_ast1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
    [kep_ast_2] = uNEO2(MJDP2,ast2,data);
    [r2,v_ast2] = sv_from_coe(kep_ast_2,ksun); % km, km/s
    [kep_ast_3] = uNEO2(MJDP3,ast3,data);
    [r3,v_ast3] = sv_from_coe(kep_ast_3,ksun); % km, km/s
    [kep_ast_4] = uNEO2(MJDP4,ast4,data);
    [r4,v_ast4] = sv_from_coe(kep_ast_4,ksun); % km, km/s % Ast 4 Arrival
    
    % Converting mJD2000 in seconds
    t1_sec = MJD01*60*60*24;
    t2a_sec = MJDP1a*60*60*24;
    t2b_sec = MJDP1b*60*60*24;
    t3_sec = MJDP2*60*60*24;
    t4_sec = MJDP3*60*60*24;
    t5_sec = MJDP4*60*60*24;
    
    %% Lamberts and deltaVs
    
    % Earth -> DSM
    [~,~,~,~,VDd,VDa,~,~] = lambertMR(rEA,rD,t1_sec,ksun,0,0,0,0);
    sol.dvD = sqrt((VDd(1)- vEA(1))^2+(VDd(2)-vEA(2))^2+(VDd(3)- vEA(3))^2);
    
    if sol.dvD < sqrt(sim.C3_max)
        sol.dV_extra_launch = 0;
        sol.Vinf_launcher = sol.dvD;
    else
        sol.dV_extra_launch = sol.dvD - sqrt(sim.C3_max);
        sol.Vinf_launcher = sqrt(sim.C3_max);
    end
    
 
    % DSM -> asteroid 1
    [~,~,~,~,V1d,V1a,~,~] = lambertMR(rD,r1,t2a_sec,ksun,0,0,0,0);
    sol.dv1 = sqrt((V1d(1)-VDa(1))^2+(V1d(2)-VDa(2))^2+(V1d(3)-VDa(3))^2);


    % asteroid 1 -> 2
    [~,~,~,~,V2d,V2a,~,~] = lambertMR(r1,r2,t2b_sec,ksun,0,0,0,0);
    sol.dv2 = sqrt((V2d(1)- V1a(1))^2+(V2d(2)- V1a(2))^2+(V2d(3)-V1a(3))^2);


    % asteroid 2 -> 3
    [~,~,~,~,V3d,V3a,~,~] = lambertMR(r2,r3,t3_sec,ksun,0,0,0,0);
    sol.dv3 = sqrt((V3d(1)-V2a(1))^2+(V3d(2)-V2a(2))^2+(V3d(3)-V2a(3))^2);


    % asteroid 3 -> 4
    [~,~,~,~,V4d,V4a,~,~] = lambertMR(r3,r4,t4_sec,ksun,0,0,0,0);
    sol.dv4 = sqrt((V4d(1)-V3a(1))^2+(V4d(2)-V3a(2))^2+(V4d(3)-V3a(3))^2);
    
    % actual total dV, different from Fval because now it enters also
    % penalty in the optimisation function
    sol.dV_tot = sol.dv1 + sol.dv2 + sol.dv3 + sol.dv4 ; 

    
    %% Plotting
    % PLOT FULL ORBITS AND BEST LAMBERT TRANSFER 
    figure('Name','Mission Orbits and Phases')
    % Earth
    plot_earth_orbit(MJD01,3,colors,8);
    hold on
    % Asteroids
    Frac_Orbit = 1/4;
    plot_asteorid_orbit(MJDP1b,Frac_Orbit,ast1,colors,2);
    plot_asteorid_orbit(MJDP2,Frac_Orbit,ast2,colors,3);
    plot_asteorid_orbit(MJDP3,Frac_Orbit,ast3,colors,4);
    plot_asteorid_orbit(MJDP4+200,Frac_Orbit,ast4,colors,5);
    
    % Mission Arcs
    % First leg: Earth -> DSM
    time_eval_12 = linspace(t1_sec,t2a_sec,sol.TOF1a);
    y012 = [rEA; VDd']; %km, km/s; velocity from lambert arc transfer orbit injection
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t12,y12] = ode113(@rates, time_eval_12, y012,options,'sun');
    ht1a = plot3( y12(:,1)./AU, y12(:,2)./AU, y12(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','Trajectory');

    % Second leg: DSM -> Ast 1
    time_eval_23 = linspace(t2a_sec,t2b_sec,sol.TOF1b);
    y023 = [rD; V1d']; %km, km/s; velocity from lambert arc transfer orbit injection
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t23,y23] = ode113(@rates, time_eval_23, y023,options,'sun');
    ht1b = plot3( y23(:,1)./AU, y23(:,2)./AU, y23(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','Trajectory');

    % Third leg: Ast 1 -> Ast 2
    time_eval_34 = linspace(t2b_sec,t3_sec,sol.TOF2);
    y034 = [r1; V2d']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t34,y34] = ode113(@rates, time_eval_34, y034,options,'sun');
    ht2 = plot3( y34(:,1)./AU, y34(:,2)./AU, y34(:,3)./AU,'Color',colors(1,:));
    ht2.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % Fourth leg: Ast 2 -> Ast 3
    time_eval_56 = linspace(t3_sec,t4_sec,sol.TOF3);
    y056 = [r2; V3d']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t56,y56] = ode113(@rates, time_eval_56, y056,options,'sun');
    ht3 = plot3( y56(:,1)./AU, y56(:,2)./AU, y56(:,3)./AU,'Color',colors(1,:));
    ht3.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % Fifth leg: Ast 3 -> Ast 4
    time_eval_78 = linspace(t4_sec,t5_sec,sol.TOF4);
    y078 = [r3; V4d']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t78,y78] = ode113(@rates, time_eval_78, y078,options,'sun');
    ht4 =plot3( y78(:,1)./AU, y78(:,2)./AU, y78(:,3)./AU,'Color',colors(1,:));
    ht4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    

    
    %% Naming
    % Sun Yellow Asterisk
    plot3(0,0,0,'*','Color',colors(4,:),'DisplayName','Sun');
    
    legend('show','Location','southeastoutside')
    
    hp1a = plot3(rEA(1)./AU,rEA(2)./AU,rEA(3)./AU,'o','Color',colors(8,:),'MarkerSize',6,...
    'DisplayName','Earth Departure');

    hp1b = plot3(rD(1)./AU,rD(2)./AU,rD(3)./AU,'o','Color',colors(10,:),'MarkerSize',6,...
    'DisplayName','DSM');
    
    hp2 = plot3(r1(1)./AU,r1(2)./AU,r1(3)./AU,'^','Color',colors(2,:),'MarkerSize',6,...
        'DisplayName',ast1);

    hp4 = plot3(r2(1)./AU,r2(2)./AU,r2(3)./AU,'^','Color',colors(3,:),'MarkerSize',6,...
        'DisplayName',ast2);

    hp6 = plot3(r3(1)./AU,r3(2)./AU,r3(3)./AU,'^','Color',colors(4,:),'MarkerSize',6,...
        'DisplayName',ast3);

    hp8 = plot3(r4(1)./AU,r4(2)./AU,r4(3)./AU,'^','Color',colors(5,:),'MarkerSize',6,...
        'DisplayName',ast4);


    axis equal; grid on
    title(sim.case_name)
    xlabel('x [AU]')
    ylabel('y [AU]')
    zlabel('z [AU]')
    
end