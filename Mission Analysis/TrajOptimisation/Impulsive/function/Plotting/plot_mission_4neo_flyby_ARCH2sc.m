function [sol] = plot_mission_4neo_flyby_ARCH2sc(sol,data,sim,colors)
    
    AU = astroConstants(2);

    %% Initialise stuff from solution
    MJD01 = sol.MJD0;
    % 1st sc
    MJDP1 = MJD01 + sol.TOF1;
    MJDP2 = MJDP1 + sol.TOF2;
    % 2nd sc
    MJDPa = MJD01 + sol.TOFa;
    MJDPb = MJDPa + sol.TOFb;
    
    ast1 = sol.ast_1;
    ast2 = sol.ast_2;
    
    ast_a = sol.ast_a;
    ast_b = sol.ast_b;
    
    %% Position of planet at solutions moments
    [kep_EA, ksun] = uplanet(MJD01, 3); % Earth Departure
    [r_EA, v_EA] = sv_from_coe(kep_EA,ksun); %km, km/s
    [kep_ast_1] = uNEO2(MJDP1,ast1,data);
    [r_ast1,v_ast1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
    [kep_ast_2] = uNEO2(MJDP2,ast2,data);
    [r_ast2,v_ast2] = sv_from_coe(kep_ast_2,ksun); % km, km/s
    
    [kep_ast_a] = uNEO2(MJDPa,ast_a,data);
    [r_asta,v_asta] = sv_from_coe(kep_ast_a,ksun); % km, km/s
    [kep_ast_b] = uNEO2(MJDPb,ast_b,data);
    [r_astb,v_astb] = sv_from_coe(kep_ast_b,ksun); % km, km/s % Ast 4 Arrival
    
    % Converting mJD2000 in seconds
    t1_sec = MJD01*60*60*24;
    
    t2_sec = MJDP1*60*60*24;
    t3_sec = MJDP2*60*60*24;
    
    ta_sec = MJDPa*60*60*24;
    tb_sec = MJDPb*60*60*24;
    
    %% Lamberts and deltaVs
    % SPACECRAFT 1
    % Earth -> asteroid 1
    ToF_EAast1_sec = t2_sec - t1_sec;
    [~,~,~,~,VI_EAast1,VF_EAast1,~,~] = lambertMR(r_EA,r_ast1,ToF_EAast1_sec,ksun,0,0,0,0);
    dv1 = sqrt((VI_EAast1(1)-v_EA(1))^2+(VI_EAast1(2)-v_EA(2))^2+(VI_EAast1(3)-v_EA(3))^2);
    if dv1 < sqrt(sim.C3_max)
        sol.dV_extra_launch_sc1 = 0;
        sol.Vinf_launcher_sc1 = dv1;
    else
        sol.dV_extra_launch_sc1 = dv1 - sqrt(sim.C3_max);
        sol.Vinf_launcher_sc1 = sqrt(sim.C3_max);
    end
    % relative velocity arrival at the asteroid and asteroid itself, for the deployment of the "lander"
    sol.Vrel_passage_ast1 = sqrt((VF_EAast1(1)-v_ast1(1))^2+(VF_EAast1(2)-v_ast1(2))^2+(VF_EAast1(3)-v_ast1(3))^2);
    
    % asteroid 1 -> asteroid 2
    ToF_ast12_sec = t3_sec - t2_sec;
    [~,~,~,~,VI_ast12,VF_ast12,~,~] = lambertMR(r_ast1,r_ast2,ToF_ast12_sec,ksun,0,0,0,0);
    % dV of flyby passage on asteroid 1
    sol.dVast1 = sqrt((VI_ast12(1)-VF_EAast1(1))^2+(VI_ast12(2)-VF_EAast1(2))^2+(VI_ast12(3)-VF_EAast1(3))^2);
    % relative velocity arrival at the asteroid and asteroid itself, for the deployment of the "lander"
    sol.Vrel_passage_ast2 = sqrt((VF_ast12(1)-v_ast2(1))^2+(VF_ast12(2)-v_ast2(2))^2+(VF_ast12(3)-v_ast2(3))^2);
    
    % SPACECRAFT 2
    % Earth -> asteroid a
    ToF_EAasta_sec = ta_sec - t1_sec;
    [~,~,~,~,VI_EAasta,VF_EAasta,~,~] = lambertMR(r_EA,r_asta,ToF_EAasta_sec,ksun,0,0,0,0);
    dva = sqrt((VI_EAasta(1)-v_EA(1))^2+(VI_EAasta(2)-v_EA(2))^2+(VI_EAasta(3)-v_EA(3))^2);
    if dva < sqrt(sim.C3_max)
        sol.dV_extra_launch_sc2 = 0;
        sol.Vinf_launcher_sc2 = dva;
    else
        sol.dV_extra_launch_sc2 = dva - sqrt(sim.C3_max);
        sol.Vinf_launcher_sc2 = sqrt(sim.C3_max);
    end
    % relative velocity arrival at the asteroid and asteroid itself, for the deployment of the "lander"
    sol.Vrel_passage_asta = sqrt((VF_EAasta(1)-v_asta(1))^2+(VF_EAasta(2)-v_asta(2))^2+(VF_EAasta(3)-v_asta(3))^2);

    % asteroid a -> asteroid b
    ToF_astab_sec = tb_sec - ta_sec;
    [~,~,~,~,VI_astab,VF_astab,~,~] = lambertMR(r_asta,r_astb,ToF_astab_sec,ksun,0,0,0,0);
    % dV of flyby passage on asteroid a
    sol.dVasta = sqrt((VI_astab(1)-VF_EAasta(1))^2+(VI_astab(2)-VF_EAasta(2))^2+(VI_astab(3)-VF_EAasta(3))^2);
    % relative velocity arrival at the asteroid and asteroid itself, for the deployment of the "lander"
    sol.Vrel_passage_astb = sqrt((VF_astab(1)-v_astb(1))^2+(VF_astab(2)-v_astb(2))^2+(VF_astab(3)-v_astb(3))^2);
    
    sol.dV_tot = sol.dV_extra_launch_sc1 + sol.dV_extra_launch_sc2 + sol.dVast1 + sol.dVasta;
    
    %% Plotting
    % PLOT FULL ORBITS AND BEST LAMBERT TRANSFER 
    figure('Name','Mission Orbits and Phases')
    % Earth
    plot_earth_orbit(MJD01,3,colors,8);
    hold on
    % Asteroids
    Frac_Orb = 1/6;
    plot_asteorid_orbit(MJDP1,Frac_Orb,ast1,colors,2);
    plot_asteorid_orbit(MJDP2+100,Frac_Orb,ast2,colors,3);
    plot_asteorid_orbit(MJDPa-50,Frac_Orb,ast_a,colors,4);
    plot_asteorid_orbit(MJDPb+50,Frac_Orb,ast_b,colors,5);
    
    % Mission Arcs
    % SC 1-First leg: Earth -> Ast 1
    t012 = t1_sec;
    tf12 = t2_sec;
    y012 = [r_EA; VI_EAast1']; %km, km/s; velocity from lambert arc transfer orbit injection
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t12,y12] = ode113(@rates, [t012 tf12], y012,options,'sun');
    plot3( y12(:,1)./AU, y12(:,2)./AU, y12(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','SC 1 Trajectory');

    % SC 1-Second leg: Ast 1 -> Ast 2
    t034 = t2_sec;
    tf34 = t3_sec;
    y034 = [r_ast1; VI_ast12']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t34,y34] = ode113(@rates, [t034 tf34], y034,options,'sun');
    hts1_2 = plot3( y34(:,1)./AU, y34(:,2)./AU, y34(:,3)./AU,'Color',colors(1,:));
    hts1_2.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % SC 2-First leg: Earth -> Ast a
    t056 = t1_sec;
    tf56 = ta_sec;
    y056 = [r_EA; VI_EAasta']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t56,y56] = ode113(@rates, [t056 tf56], y056,options,'sun');
    plot3( y56(:,1)./AU, y56(:,2)./AU, y56(:,3)./AU,'Color',colors(10,:),...
        'DisplayName','SC 2 Trajectory');

    % SC 2-Fourth leg: Ast a -> Ast b
    t078 = ta_sec;
    tf78 = tb_sec;
    y078 = [r_asta; VI_astab']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t78,y78] = ode113(@rates, [t078 tf78], y078,options,'sun');
    hts2_2 = plot3( y78(:,1)./AU, y78(:,2)./AU, y78(:,3)./AU,'Color',colors(10,:));
    hts2_2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    %% Extract the Sun-SpaceCraft Distance for all the trajectory [km]
    sol.SunSpacecraftDistanceNorm = vecnorm([y12;y34;y56;y78],2,2); % 2,2 means norm 2 and by row
    sol.SpacecraftTrajectory = [y12;y34;y56;y78];
    sol.SCtime = [t12;t34;t56;t78];
    [sol.angles.SAA,sol.angles.EVA,sol.angles.SCA,sol.angles.SolarConjunction] = aspect_angles(sol);
    
    %% Naming
    % Sun Yellow Asterisk
    plot3(0,0,0,'*','Color',colors(4,:),'DisplayName','Sun');
    
    legend('show','Location','southeastoutside')
    
    hp1 = plot3(r_EA(1)./AU,r_EA(2)./AU,r_EA(3)./AU,'o','Color',colors(8,:),'MarkerSize',6,...
    'DisplayName','Earth Departure');
%     hp1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp2 = plot3(r_ast1(1)./AU,r_ast1(2)./AU,r_ast1(3)./AU,'^','Color',colors(2,:),'MarkerSize',6,...
        'DisplayName',ast1);
%     hp2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp4 = plot3(r_ast2(1)./AU,r_ast2(2)./AU,r_ast2(3)./AU,'^','Color',colors(3,:),'MarkerSize',6,...
        'DisplayName',ast2);
%     hp4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp6 = plot3(r_asta(1)./AU,r_asta(2)./AU,r_asta(3)./AU,'^','Color',colors(4,:),'MarkerSize',6,...
        'DisplayName',ast_a);
%     hp6.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp8 = plot3(r_astb(1)./AU,r_astb(2)./AU,r_astb(3)./AU,'^','Color',colors(5,:),'MarkerSize',6,...
        'DisplayName',ast_b);
%     hp8.Annotation.LegendInformation.IconDisplayStyle = 'off';

    axis equal; grid on
%     title(sim.case_name)
    xlabel('x [AU]')
    ylabel('y [AU]')
    zlabel('z [AU]')
    
end