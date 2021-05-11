function [sol] = plot_mission_4neo_RV_ARCH2sc(sol,data,sim,colors)
    
    AU = astroConstants(2);

    %% Initialise stuff from solution	
    MJD01 = sol.MJD0;
    % 1st sc
    MJDF1 = MJD01 + sol.TOF1;
    MJD02 = MJDF1 + sol.buffer_time1; % departure of sc 1 from ast 1
    MJDF2 = MJD02 + sol.TOF2;
    
    % 2nd sc
    MJDFa = MJD01 + sol.TOFa;
    MJD0a = MJDFa + sol.buffer_time2; % departure of sc 2 from ast a
    MJDFb = MJD0a + sol.TOFb;
    
    ast1 = sol.ast_1;
    ast2 = sol.ast_2;
    
    ast_a = sol.ast_a;
    ast_b = sol.ast_b;
    
    %% Position of planet at solutions moments
    [kep_EA, ksun] = uplanet(MJD01, 3); % Earth Departure
    [r_EA, v_EA] = sv_from_coe(kep_EA,ksun); %km, km/s
    % arrival of 1st sc at 1st ast
    [kep_ast_11] = uNEO2(MJDF1,ast1,data); % [km,-,rad,rad,rad,wrapped rad]
    [r1_arr, v1_arr] = sv_from_coe(kep_ast_11,ksun); % km, km/s
    % departure of 1st sc at 1st ast
    [kep_ast_12] = uNEO2(MJD02,ast1,data); % [km,-,rad,rad,rad,wrapped rad]
    [r1_dep, v1_dep] = sv_from_coe(kep_ast_12,ksun); % km, km/s
    [kep_ast_2] = uNEO2(MJDF2,ast2,data);
    [r_ast2,v_ast2] = sv_from_coe(kep_ast_2,ksun); % km, km/s
    
    % arrival of 2nd sc at 1st ast
    [kep_ast_a1] = uNEO2(MJDFa,ast_a,data);
    [ra_arr, va_arr] = sv_from_coe(kep_ast_a1,ksun); % km, km/s
    % departure of 2nd sc at 1st ast
    [kep_ast_a2] = uNEO2(MJD0a,ast_a,data);
    [ra_dep, va_dep] = sv_from_coe(kep_ast_a2,ksun); % km, km/s
    [kep_ast_b] = uNEO2(MJDFb,ast_b,data);
    [r_astb,v_astb] = sv_from_coe(kep_ast_b,ksun); % km, km/s % Ast b Arrival
    
    % Converting mJD2000 in seconds
    ToF_EAast1_sec = sol.TOF1*60*60*24;
    buffer_time1_sec = sol.buffer_time1*60*60*24;
    ToF_ast12_sec = sol.TOF2*60*60*24;

    ToF_EAasta_sec = sol.TOFa*60*60*24;
    buffer_time2_sec = sol.buffer_time2*60*60*24;
    ToF_astab_sec = sol.TOFb*60*60*24;
    
    %% Lamberts and deltaVs
    % SPACECRAFT 1
    % Earth -> asteroid 1
    [~,~,~,~,VI_EAast1,VF_EAast1,~,~] = lambertMR(r_EA,r1_arr,ToF_EAast1_sec,ksun,0,0,0,0);
    dv1 = sqrt((VI_EAast1(1)-v_EA(1))^2+(VI_EAast1(2)-v_EA(2))^2+(VI_EAast1(3)-v_EA(3))^2);
    if dv1 < sqrt(sim.C3_max)
        sol.dV_extra_launch_sc1 = 0;
        sol.Vinf_launcher_sc1 = dv1;
    else
        sol.dV_extra_launch_sc1 = dv1 - sqrt(sim.C3_max);
        sol.Vinf_launcher_sc1 = sqrt(sim.C3_max);
    end
    
    sol.dV2_EAast1_sc1 = sqrt((VF_EAast1(1)-v1_arr(1))^2+(VF_EAast1(2)-v1_arr(2))^2+(VF_EAast1(3)-v1_arr(3))^2);
%     sol.dV_tot_sc1_leg1 = sol.dV_EAast1_sc1 + sol.dV_extra_launch_sc1;
    
    % asteroid 1 -> asteroid 2
    [~,~,~,~,VI_ast12,VF_ast12,~,~] = lambertMR(r1_dep,r_ast2,ToF_ast12_sec,ksun,0,0,0,0);
    sol.dV1_ast12 = sqrt((VI_ast12(1)-v1_dep(1))^2+(VI_ast12(2)-v1_dep(2))^2+(VI_ast12(3)-v1_dep(3))^2);
    sol.dV2_ast12 = sqrt((VF_ast12(1)-v_ast2(1))^2+(VF_ast12(2)-v_ast2(2))^2+(VF_ast12(3)-v_ast2(3))^2);
    % dV total of rendezvous on asteroid 1, enter the asteroid orbit and exit
    % to go toward asteroid 2
    sol.dV_RV_sc1 = sol.dV_extra_launch_sc1 + sol.dV2_EAast1_sc1 + sol.dV1_ast12 + sol.dV2_ast12;

    % SPACECRAFT 2
    % Earth -> asteroid a
    [~,~,~,~,VI_EAasta,VF_EAasta,~,~] = lambertMR(r_EA,ra_arr,ToF_EAasta_sec,ksun,0,0,0,0);
    dva = sqrt((VI_EAasta(1)-v_EA(1))^2+(VI_EAasta(2)-v_EA(2))^2+(VI_EAasta(3)-v_EA(3))^2);
    if dva < sqrt(sim.C3_max)
        sol.dV_extra_launch_sc2 = 0;
        sol.Vinf_launcher_sc2 = dva;
    else
        sol.dV_extra_launch_sc2 = dva - sqrt(sim.C3_max);
        sol.Vinf_launcher_sc2 = sqrt(sim.C3_max);
    end
    sol.dV2_EAasta_sc2 = sqrt((VF_EAasta(1)-va_arr(1))^2+(VF_EAasta(2)-va_arr(2))^2+(VF_EAasta(3)-va_arr(3))^2);
    
    % asteroid a -> asteroid b
    [~,~,~,~,VI_astab,VF_astab,~,~] = lambertMR(ra_dep,r_astb,ToF_astab_sec,ksun,0,0,0,0);
     % asteroid 1 -> asteroid 2
    sol.dV1_astab = sqrt((VI_astab(1)-va_dep(1))^2+(VI_astab(2)-va_dep(2))^2+(VI_astab(3)-va_dep(3))^2);
    sol.dV2_astab = sqrt((VF_astab(1)-v_astb(1))^2+(VF_astab(2)-v_astb(2))^2+(VF_astab(3)-v_astb(3))^2);
    
    % dV total of rendezvous on asteroid a, enter the asteroid orbit and exit
    % to go toward asteroid b
    sol.dV_RV_sc2 = sol.dV_extra_launch_sc2 + sol.dV2_EAasta_sc2 + sol.dV1_astab + sol.dV2_astab;
    sol.dV_tot = sol.dV_RV_sc1 + sol.dV_RV_sc2;
    
    %% Plotting
    % PLOT FULL ORBITS AND BEST LAMBERT TRANSFER 
    figure('Name','Mission Orbits and Phases')
    % Earth
    plot_earth_orbit(MJD01,3,colors,8);
    hold on
    % Asteroids
    frac_orbit = 1/6;
    plot_asteorid_orbit(MJDF1,frac_orbit,ast1,colors,2); % ast1
    plot_asteorid_orbit(MJDF2-500,frac_orbit,ast2,colors,3); % ast2
    
    plot_asteorid_orbit(MJDFa,frac_orbit,ast_a,colors,4); % asta
    plot_asteorid_orbit(MJDFb,frac_orbit,ast_b,colors,5); % astb
    
    %% Mission Arcs
    % SC 1-First leg: Earth -> Ast 1
    departure_earth_sec = MJD01*3600*24;
    arrival_ast1_sec = MJDF1*3600*24;
    y0_EAast1 = [r_EA; VI_EAast1']; %km, km/s; velocity from lambert arc transfer orbit injection
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [tEA1,yEA1] = ode113(@rates, [departure_earth_sec arrival_ast1_sec], y0_EAast1,options,'sun');
    plot3( yEA1(:,1)./AU, yEA1(:,2)./AU, yEA1(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','SC 1 Trajectory');

    % Coasting on Ast 1
    start_ast1_sec = arrival_ast1_sec;
    end_ast1_sec = MJD02*3600*24;
%     time_int23 = [start_ast1_sec, end_ast1_sec];
    time_int23 = linspace(start_ast1_sec,end_ast1_sec,50);
    y0_coast_ast1 = [r1_arr; v1_arr]; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [tC1,yC1] = ode113(@rates, time_int23,y0_coast_ast1,options,'sun');
    hc1 = plot3( yC1(:,1)./AU, yC1(:,2)./AU, yC1(:,3)./AU,'*','Color',colors(1,:),'Markersize',3);
    hc1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    % SC 1-Second leg: Ast 1 -> Ast 2
    departure_ast1_sec = end_ast1_sec;
    arrival_ast2_sec = MJDF2*3600*24;
    y0_ast12 = [r1_dep; VI_ast12']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t_ast12,y_ast12] = ode113(@rates, [departure_ast1_sec arrival_ast2_sec], y0_ast12,options,'sun');
    hts1_2 = plot3( y_ast12(:,1)./AU, y_ast12(:,2)./AU, y_ast12(:,3)./AU,'Color',colors(1,:));
    hts1_2.Annotation.LegendInformation.IconDisplayStyle = 'off';

    %% SC 2-First leg: Earth -> Ast a
    arrival_asta_sec = MJDFa*3600*24;
    y0_EA_asta = [r_EA; VI_EAasta']; %km, km/s
    time_int_EA_asta = [departure_earth_sec, arrival_asta_sec];
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t_EAa,y_EAa] = ode113(@rates, time_int_EA_asta, y0_EA_asta,options,'sun');
    plot3( y_EAa(:,1)./AU, y_EAa(:,2)./AU, y_EAa(:,3)./AU,'Color',colors(10,:),...
        'DisplayName','SC 2 Trajectory');
    
    % Coasting on Ast a
    start_asta_sec = arrival_asta_sec;
    end_asta_sec = MJD0a*3600*24;
    time_int_Ca = linspace(start_asta_sec,end_asta_sec,50);
    y0Ca = [ra_arr; va_arr]; %km, km/s; after the application of 2nd dv of its transfer leg, rendezvous with asteroid, same velocity as the asteroid at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [tCa,yCa] = ode113(@rates, time_int_Ca, y0Ca,options,'sun');
    hc2 = plot3( yCa(:,1)./AU, yCa(:,2)./AU, yCa(:,3)./AU,'*','Color',colors(10,:),'Markersize',3);
    hc2.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % SC 2-Fourth leg: Ast a -> Ast b
    departure_asta_sec = end_asta_sec;
    arrival_astb_sec = MJDFb*3600*24;
    time_int_astab = [departure_asta_sec,arrival_astb_sec];
    y0_astab = [ra_dep; VI_astab']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t_astab,y_astab] = ode113(@rates, time_int_astab, y0_astab,options,'sun');
    hts2_2 = plot3( y_astab(:,1)./AU, y_astab(:,2)./AU, y_astab(:,3)./AU,'Color',colors(10,:));
    hts2_2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    %% Extract the Sun-SpaceCraft Distance for all the trajectory [km]
%     sol.SunSpacecraftDistanceNorm = vecnorm([yEA1;y_ast12;y_EAa;y_astab],2,2); % 2,2 means norm 2 and by row
    sol.Spacecraft1Trajectory = [yEA1;yC1;y_ast12];
    sol.Spacecraft2Trajectory = [y_EAa;yCa;y_astab];
    sol.SC1time = [tEA1;tC1;t_ast12;];
    sol.SC2time = [t_EAa;tCa;t_astab];
%     [sol.angles.SAA,sol.angles.EVA,sol.angles.SCA,sol.angles.SolarConjunction] = aspect_angles(sol);
    
    %% Naming
    % Sun Yellow Asterisk
    plot3(0,0,0,'*','Color',colors(4,:),'DisplayName','Sun');
    
    legend('show','Location','southeastoutside')
    
    hp1 = plot3(r_EA(1)./AU,r_EA(2)./AU,r_EA(3)./AU,'o','Color',colors(8,:),'MarkerSize',6,...
    'DisplayName','Earth Departure');
%     hp1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp2 = sprintf('Arrival %s',ast1);
    hp2 = plot3(r1_arr(1)./AU,r1_arr(2)./AU,r1_arr(3)./AU,'^','Color',colors(2,:),'MarkerSize',6,...
    'DisplayName',name_hp2);
    name_hp3 = sprintf('Departure %s',ast1);
    hp3 = plot3(r1_dep(1)./AU,r1_dep(2)./AU,r1_dep(3)./AU,'o','Color',colors(2,:),'MarkerSize',6,...
    'DisplayName',name_hp3);
    name_hp4 = sprintf('Arrival %s',ast2);
    hp4 = plot3(r_ast2(1)./AU,r_ast2(2)./AU,r_ast2(3)./AU,'^','Color',colors(3,:),'MarkerSize',6,...
    'DisplayName',name_hp4);

    name_hp5 = sprintf('Arrival %s',ast_a);
    hp5 = plot3(ra_arr(1)./AU,ra_arr(2)./AU,ra_arr(3)./AU,'o','Color',colors(4,:),'MarkerSize',6,...
    'DisplayName',name_hp5);
    name_hp6 = sprintf('Departure %s',ast_a);
    hp6 = plot3(ra_dep(1)./AU,ra_dep(2)./AU,ra_dep(3)./AU,'^','Color',colors(4,:),'MarkerSize',6,...
    'DisplayName',name_hp6);
%     hp6.Annotation.LegendInformation.IconDisplayStyle = 'off';
    name_hp7 = sprintf('Arrival %s',ast_b);
    hp7 = plot3(r_astb(1)./AU,r_astb(2)./AU,r_astb(3)./AU,'o','Color',colors(5,:),'MarkerSize',6,...
    'DisplayName',name_hp7);

    axis equal; grid on
    title(sim.case_name)
    xlabel('x [AU]')
    ylabel('y [AU]')
    zlabel('z [AU]')
    
end