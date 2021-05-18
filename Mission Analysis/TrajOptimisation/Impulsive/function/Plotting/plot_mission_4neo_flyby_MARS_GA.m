function [sol] = plot_mission_4neo_flyby_MARS_GA(sol,asteroid_names_sequence,data,sim,colors)
    
    AU = astroConstants(2);
    
    %% Initialise stuff from solution
    MJD01 = sol.MJD0;
    MJDGA = MJD01 + sol.TOF0;
    MJDP1 = MJDGA + sol.TOF1;
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
    [kep_GA, ksun] = uplanet(MJDGA,4); % Mars Gravity Assist
    [rGA, vGA] = sv_from_coe(kep_GA,ksun); %km, km/s
    [kep_ast_1] = uNEO2(MJDP1,ast1,data);
    [r1,v1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
    [kep_ast_2] = uNEO2(MJDP2,ast2,data);
    [r2,v2] = sv_from_coe(kep_ast_2,ksun); % km, km/s
    [kep_ast_3] = uNEO2(MJDP3,ast3,data);
    [r3,v3] = sv_from_coe(kep_ast_3,ksun); % km, km/s
    [kep_ast_4] = uNEO2(MJDP4,ast4,data);
    [r4,v4] = sv_from_coe(kep_ast_4,ksun); % km, km/s % Ast 4 Arrival
    
    % Converting mJD2000 in seconds
    ToF_EAGA_sec = sol.TOF0*60*60*24;
    ToF_GAast1_sec = sol.TOF1*60*60*24;
    ToF_ast12_sec = sol.TOF2*60*60*24;
    ToF_ast23_sec = sol.TOF3*60*60*24;
    ToF_ast34_sec = sol.TOF4*60*60*24;
    
    %% Lamberts and deltaVs
    [~,~,~,~,VI_EAGA,VF_EAGA,~,~] = lambertMR(rEA,rGA,ToF_EAGA_sec,ksun,0,0,0,0);
    dv1 = sqrt((VI_EAGA(1)-vEA(1))^2+(VI_EAGA(2)-vEA(2))^2+(VI_EAGA(3)-vEA(3))^2);
    if dv1 < sqrt(sim.C3_max)
        sol.dV_extra_launch = 0;
        sol.Vinf_launcher = dv1;
    else
        sol.dV_extra_launch = dv1 - sqrt(sim.C3_max);
        sol.Vinf_launcher = sqrt(sim.C3_max);
    end
%     sol.dV_tot_leg1 = sol.dV_single.dV2 + sol.dV_single.dV_extra_launch;
    % relative velocity arrival at the asteroid and asteroid itself, for the deployment of the "lander"
%     sol.Vrel_passage_ast1 = sqrt((VF_EAGA(1)-v_ast1(1))^2+(VF_EAGA(2)-v_ast1(2))^2+(VF_EAGA(3)-v_ast1(3))^2);
    [~,~,~,~,VI_GAast1,VF_GAast1,~,~] = lambertMR(rGA,r1,ToF_GAast1_sec,ksun,0,0,0,0);
    dv2_GAast1 = sqrt((VF_GAast1(1)-v1(1))^2+(VF_GAast1(2)-v1(2))^2+(VF_GAast1(3)-v1(3))^2);
    
    if sim.ID_FLYBY == 3
        RPlanet_flyby = astroConstants(23); % Radius_Earth, km
        muPlanet_flyby = astroConstants(13); % muEarth, km^3/s^2
        R_lim_from_planet = 500; % km, for earth is ok to avoid atmosphere
    elseif sim.ID_FLYBY == 4
        RPlanet_flyby = astroConstants(24); % Radius_mars, km
        muPlanet_flyby = astroConstants(14); % mu mars, km^3/s^2
        R_lim_from_planet = 200; % km, for mars is ok to avoid atmosphere
    end
     sol.delta_V_p = flyby(RPlanet_flyby, muPlanet_flyby,R_lim_from_planet, ...
                  MJDGA, VF_EAGA', VI_GAast1', sim.ID_FLYBY); % input needs to be vertical vectors

    % Velocity gained with flyby
    sol.dV_gained_flyby = sqrt((VI_GAast1(1)-VF_EAGA(1))^2+(VI_GAast1(2)-VF_EAGA(2))^2+(VI_GAast1(3)-VF_EAGA(3))^2) - sol.delta_V_p;
    sol.Vrel_passage_ast1 = sqrt((VF_GAast1(1)-v1(1))^2+(VF_GAast1(2)-v1(2))^2+(VF_GAast1(3)-v1(3))^2);
    
    % asteroid 1 -> 2
    [~,~,~,~,VI_ast12,VF_ast12,~,~] = lambertMR(r1,r2,ToF_ast12_sec,ksun,0,0,0,0);
    dv2_ast12 = sqrt((VF_ast12(1)-v2(1))^2+(VF_ast12(2)-v2(2))^2+(VF_ast12(3)-v2(3))^2);
    
    % dV of flyby passage on asteroid 1
    sol.dVast1 = sqrt((VI_ast12(1)-VF_GAast1(1))^2+(VI_ast12(2)-VF_GAast1(2))^2+(VI_ast12(3)-VF_GAast1(3))^2);
    % relative velocity arrival at the asteroid and asteroid itself, for the deployment of the "lander"
    sol.Vrel_passage_ast2 = sqrt((VF_ast12(1)-v2(1))^2+(VF_ast12(2)-v2(2))^2+(VF_ast12(3)-v2(3))^2);
    
    % asteroid 2 -> 3
    [~,~,~,~,VI_ast23,VF_ast23,~,~] = lambertMR(r2,r3,ToF_ast23_sec,ksun,0,0,0,0);
    dv2_ast23 = sqrt((VF_ast23(1)-v3(1))^2+(VF_ast23(2)-v3(2))^2+(VF_ast23(3)-v3(3))^2);

    % dV of flyby passage on asteroid 2
    sol.dVast2 = sqrt((VI_ast23(1)-VF_ast12(1))^2+(VI_ast23(2)-VF_ast12(2))^2+(VI_ast23(3)-VF_ast12(3))^2);
    % relative velocity arrival at the asteroid and asteroid itself, for the deployment of the "lander"
    sol.Vrel_passage_ast3 = sqrt((VF_ast23(1)-v3(1))^2+(VF_ast23(2)-v3(2))^2+(VF_ast23(3)-v3(3))^2);
    
    % asteorid 3 -> 4
    [~,~,~,~,VI_ast34,VF_ast34,~,~] = lambertMR(r3,r4,ToF_ast34_sec,ksun,0,0,0,0);
    dv2_ast34 = sqrt((VF_ast34(1)-v4(1))^2+(VF_ast34(2)-v4(2))^2+(VF_ast34(3)-v4(3))^2);

    % dV of flyby passage on asteroid 3
    sol.dVast3 = sqrt((VI_ast34(1)-VF_ast23(1))^2+(VI_ast34(2)-VF_ast23(2))^2+(VI_ast34(3)-VF_ast23(3))^2);
    % relative velocity arrival at the asteroid and asteroid itself, for the deployment of the "lander"
    sol.Vrel_passage_ast4 = sqrt((VF_ast34(1)-v4(1))^2+(VF_ast34(2)-v4(2))^2+(VF_ast34(3)-v4(3))^2);
    
    % actual total dV, different from Fval because now it enters also
    % penalty in the optimisation function
    sol.dV_tot = sol.dV_extra_launch + sol.dVast1 + sol.dVast2 + sol.dVast3 + sol.delta_V_p;
    
    %% Plots
    % PLOT FULL ORBITS AND BEST LAMBERT TRANSFER 
    figure('Name','Mission Orbits and Phases')
    % Planets
    plot_earth_orbit(MJD01,3,colors,8); % Earth
    hold on
    plot_earth_orbit(MJDGA,4,colors,6); % Mars
    % Asteroids
    Orbit_Fraction = 1/6;
    plot_asteorid_orbit(MJDP1,Orbit_Fraction,ast1,colors,2);
    plot_asteorid_orbit(MJDP2,Orbit_Fraction,ast2,colors,3);
    plot_asteorid_orbit(MJDP3+200,Orbit_Fraction,ast3,colors,4);
    plot_asteorid_orbit(MJDP4,Orbit_Fraction,ast4,colors,5);
    
    % Mission Arcs
    % Converting mJD2000 passage time in seconds
    tEA_sec = MJD01*60*60*24;
    tGA_sec = MJDGA*60*60*24;
    t1_sec = MJDP1*60*60*24;
    t2_sec = MJDP2*60*60*24;
    t3_sec = MJDP3*60*60*24;
    t4_sec = MJDP4*60*60*24;
    
    % First leg: Earth -> Mars GA
    y0EAGA = [rEA; VI_EAGA']; %km, km/s; velocity from lambert arc transfer orbit injection
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [tEAGA,yEAGA] = ode113(@rates, [tEA_sec tGA_sec], y0EAGA,options,'sun');
    plot3( yEAGA(:,1)./AU, yEAGA(:,2)./AU, yEAGA(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','Trajectory');
    
    % second leg: Mars -> Ast 1
    y0GAast1 = [rGA; VI_GAast1']; %km, km/s; velocity from lambert arc transfer orbit injection
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [tEAast1,yGAast1] = ode113(@rates, [tGA_sec t1_sec], y0GAast1,options,'sun');
    ht1 = plot3( yGAast1(:,1)./AU, yGAast1(:,2)./AU, yGAast1(:,3)./AU,'Color',colors(1,:));
    ht1.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % 3rd leg: Ast 1 -> Ast 2
    y034 = [r1; VI_ast12']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t34,y34] = ode113(@rates, [t1_sec t2_sec], y034,options,'sun');
    ht2 = plot3( y34(:,1)./AU, y34(:,2)./AU, y34(:,3)./AU,'Color',colors(1,:));
    ht2.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % 4th leg: Ast 2 -> Ast 3
    y056 = [r2; VI_ast23']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t56,y56] = ode113(@rates, [t2_sec t3_sec], y056,options,'sun');
    ht3 = plot3( y56(:,1)./AU, y56(:,2)./AU, y56(:,3)./AU,'Color',colors(1,:));
    ht3.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % 5th leg: Ast 3 -> Ast 4
    y078 = [r3; VI_ast34']; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [t78,y78] = ode113(@rates, [t3_sec t4_sec], y078,options,'sun');
    ht4 = plot3( y78(:,1)./AU, y78(:,2)./AU, y78(:,3)./AU,'Color',colors(1,:));
    ht4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    %% Extract the Sun-SpaceCraft Distance for all the trajectory [km]
    sol.SunSpacecraftDistanceNorm = vecnorm([yEAGA;yGAast1;y34;y56;y78],2,2); % 2,2 means norm 2 and by row
    sol.SpacecraftTrajectory = [yEAGA;yGAast1;y34;y56;y78];
    sol.SCtime = [tEAGA;tEAast1;t34;t56;t78];
    [sol.angles.SAA,sol.angles.EVA,sol.angles.SCA,sol.angles.SolarConjunction] = aspect_angles(sol);
    
    %% Naming
    % Sun Yellow Asterisk
    plot3(0,0,0,'*','Color',colors(4,:),'DisplayName','Sun');
    
    legend('show','Location','southeastoutside')
    
    hp1 = plot3(rEA(1)./AU,rEA(2)./AU,rEA(3)./AU,'o','Color',colors(8,:),'MarkerSize',6,...
        'DisplayName','Earth Departure');
%     hp1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp1 = plot3(rGA(1)./AU,rGA(2)./AU,rGA(3)./AU,'^','Color',colors(6,:),'MarkerSize',6,...
        'DisplayName','Mars GA');
%     hp1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp2 = plot3(r1(1)./AU,r1(2)./AU,r1(3)./AU,'^','Color',colors(2,:),'MarkerSize',6,...
        'DisplayName',ast1);
%     hp2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp4 = plot3(r2(1)./AU,r2(2)./AU,r2(3)./AU,'^','Color',colors(3,:),'MarkerSize',6,...
        'DisplayName',ast2);
%     hp4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp6 = plot3(r3(1)./AU,r3(2)./AU,r3(3)./AU,'^','Color',colors(4,:),'MarkerSize',6,...
        'DisplayName',ast3);
%     hp6.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp8 = plot3(r4(1)./AU,r4(2)./AU,r4(3)./AU,'^','Color',colors(5,:),'MarkerSize',6,...
        'DisplayName',ast4);
%     hp8.Annotation.LegendInformation.IconDisplayStyle = 'off';

    axis equal; grid on
%     title(sim.case_name)
    xlabel('x [AU]')
    ylabel('y [AU]')
    zlabel('z [AU]')
    
    %% Mass Calculations
    g0 = sim.g0; %km/s^2
    Isp_mother=sim.Isp_mother; %s
    Isp_lander=sim.Isp_lander; %s
    dry_mass=sim.dry_mass; %kg
    dry_mass_lander=sim.dry_mass_lander; %kg
    % m_wet/m_dry = exp(dv/(Isp*g0))
    % mass of the lander required to pay to cancel dVrel at each flyby passage
    % dry mass lander and the propellant related to that dV
    mass_lander=[dry_mass_lander*exp(dv2_GAast1/(Isp_lander*g0)), dry_mass_lander*exp(dv2_ast12/(Isp_lander*g0)),...
                 dry_mass_lander*exp(dv2_ast23/(Isp_lander*g0)), dry_mass_lander*exp(dv2_ast34/(Isp_lander*g0))];
    % Back interpolate the mass of the overall spacecraft from end of mission,
    % with an end mass of sim.dry_mass and we build back the initial wet mass
    % at each asteroid encounter, the mothercraft expell a lander
    sol.mass.final_dry_mass = dry_mass;
    sol.mass.bef_ast4=dry_mass+mass_lander(4);
    % at each flyby the sc expell mass to perform the impulsive dV
    sol.mass.aft_ast3=sol.mass.bef_ast4*exp((sol.dVast3)/(Isp_mother*g0));
    sol.mass.bef_ast3=sol.mass.aft_ast3+mass_lander(3);
    sol.mass.aft_ast2=sol.mass.bef_ast3*exp((sol.dVast2)/(Isp_mother*g0));
    sol.mass.bef_ast2=sol.mass.aft_ast2+mass_lander(2);
    sol.mass.aft_ast1=sol.mass.bef_ast2*exp((sol.dVast1)/(Isp_mother*g0));
    sol.mass.bef_ast1=sol.mass.aft_ast1+mass_lander(1);
    % mass expelled from powered gravity assist on the earth
    sol.mass.aft_ga = sol.mass.bef_ast1*exp((sol.delta_V_p)/(Isp_mother*g0));
    % mass used due to the assisted insertion in the 1st interplanetary leg,
    % other then the sqrt(C3) of the launcher, if any
    sol.mass.departure_wet_mass=sol.mass.aft_ga*exp((sol.dV_extra_launch)/(Isp_mother*g0));

end