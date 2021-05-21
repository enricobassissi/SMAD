function obj_fun = ff_ea_4ast_LT_GA_FB_soo_NLI_ifless(x,sim,data)
%{ 
input variable vector
x = [(1) MJD0,
     (2) TOF1,
     (3) TOF2,
     (4) TOF3,
     (5) TOF4
     (6) NREV,
     (7) NREV2,
     (8) NREV3,
     (9) NREV4
     (10) IDP,
     (11) v_inf_magn,
     (12) azimuth,
     (13) elevation,
     (14) TOFGA
     GA Stuff in
     v_inf_magn2 (15)
     azimuth2 (16)
     elevation2 (17)
     GA Stuff out
     azimuth3 (18)
     elevation3 (19)
	 N REV GA (20)]
%}

%% setting the input times
MJD01 = x(1); % departure time from earth
TOFGA = x(14);
MJDPGA = MJD01 + TOFGA; % passage time GA
TOF1 = x(2);
MJDP1 = MJD01 + TOF1; % passage time on ast 1
TOF2 = x(3);
MJDP2 = MJDP1 + TOF2; % passage ast 2
TOF3 = x(4);
MJDP3 = MJDP2 + TOF3; % passage ast 3
TOF4 = x(5);
MJDP4 = MJDP3 + TOF4; % passage ast 4

    
    % N REV1
    N_rev1 = x(6);
    % N REV2
    N_rev2 = x(7);
    % N REV3
    N_rev3 = x(8);
    % N REV4
    N_rev4 = x(9);
    % N REV GA
    N_revGA = x(20);

    % C3 launcher
    v_inf_magn = x(11);
    az = x(12);
    elev = x(13);
    
    % GA stuff IN
    v_inf_magn2 = x(15);
    az2 = x(16);
    el2 = x(17);
    % GA stuff out
    az3 = x(18);
    el3 = x(19);

    % chosing which asteroid to visit
    IDP = x(10); %index of permutation, the column of the Permutation Matrix of the asteroids
    asteroid_1 = data.PermutationMatrix(IDP,1);
    asteroid_2 = data.PermutationMatrix(IDP,2);
    asteroid_3 = data.PermutationMatrix(IDP,3);
    asteroid_4 = data.PermutationMatrix(IDP,4);

    %% Computing position and velocity of the planets in that days
    % Departure from Earth
    MJD01_dim = MJD01*sim.TU/(3600*24);
    [kep_EA,ksun] = uplanet(MJD01_dim, 3);
    [r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
    r_EA = r_EA/sim.DU;
    v_EA = v_EA/sim.DU*sim.TU;

    % GA
    MJDPGA_dim = MJDPGA*sim.TU/(3600*24);
    [kep_GA,ksun] = uplanet(MJDPGA_dim, sim.ID_FLYBY);
    [r_GA, v_GA] = sv_from_coe(kep_GA,ksun);
    r_GA = r_GA/sim.DU;
    v_GA = v_GA/sim.DU*sim.TU;
    
    % PASSAGE at 1st ast
    MJDP1_dim = MJDP1*sim.TU/(3600*24);
    [kep_ast_1] = uNEO2(MJDP1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
    [r1, v1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
    r1 = r1/sim.DU;
    v1 = v1/sim.DU*sim.TU;

    % passage at 2nd ast
    MJDP2_dim = MJDP2*sim.TU/(3600*24);
    [kep_ast_2] = uNEO2(MJDP2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
    [r2, v2] = sv_from_coe(kep_ast_2,ksun); % km, km/s
    r2 = r2/sim.DU;
    v2 = v2/sim.DU*sim.TU;

    % passage at 3rd ast
    MJDP3_dim = MJDP3*sim.TU/(3600*24);
    [kep_ast_3] = uNEO2(MJDP3_dim,asteroid_3,data); % [km,-,rad,rad,rad,wrapped rad]
    [r3, v3] = sv_from_coe(kep_ast_3,ksun); % km, km/s
    r3 = r3/sim.DU;
    v3 = v3/sim.DU*sim.TU;

    % passage at 4th ast
    MJDP4_dim = MJDP4*sim.TU/(3600*24);
    [kep_ast_4] = uNEO2(MJDP4_dim,asteroid_4,data); % [km,-,rad,rad,rad,wrapped rad]
    [r4, v4] = sv_from_coe(kep_ast_4,ksun); % km, km/s
    r4 = r4/sim.DU;
    v4 = v4/sim.DU*sim.TU;

    %% Launcher departure variable
%     z = r .* sin(elev);
%     rcoselev = r .* cos(elev);
%     x = rcoselev .* cos(az);
%     y = rcoselev .* sin(az);
    v_launcher = v_inf_magn*[cos(elev)*cos(az); cos(elev)*sin(az); sin(elev)];
    v_dep = v_EA + v_launcher;  %if parabolic escape (v_extra = 0)
    
    % Gravity assist
    v_arr_GA = v_GA + v_inf_magn2*[cos(el2)*cos(az2); cos(el2)*sin(az2); sin(el2)];
   
    %% NLI
    obj_fun = 10;
    % 1st leg - Earth -> GA
    [output_GA] = NL_interpolator( r_EA , r_GA , v_dep , v_arr_GA , N_revGA , TOFGA , sim.M ,sim.PS.Isp , sim );
    % 2nd leg - GA -> Ast 1
    v_dep_GA = v_GA + v_inf_magn2*[cos(el3)*cos(az3); cos(el3)*sin(az3); sin(el3)];
    if sim.ID_FLYBY == 3
        RPlanet_flyby = astroConstants(23); % Radius_Earth, km
        muPlanet_flyby = astroConstants(13); % muEarth, km^3/s^2
        R_lim_from_planet = 500; % km, for earth is ok to avoid atmosphere
        R_SOI_PL = 0.929*1e6; % km
    elseif sim.ID_FLYBY == 4
        RPlanet_flyby = astroConstants(24); % Radius_Mars, km
        muPlanet_flyby = astroConstants(14); % mumarte, km^3/s^2
        R_lim_from_planet = 200; % km, for mars is ok to avoid atmosphere
        R_SOI_PL = 0.578*1e6; %km
    end
    MJDPGA_dim = MJDPGA*sim.TU/(3600*24);
    v_arr_GA_dim = v_arr_GA.*sim.DU./sim.TU;
    v_dep_GA_dim = v_dep_GA.*sim.DU./sim.TU;
    [delta_v_p, rp] = flyby(RPlanet_flyby, muPlanet_flyby,R_lim_from_planet, ...
                      MJDPGA_dim, v_arr_GA_dim, v_dep_GA_dim, sim.ID_FLYBY, R_SOI_PL);
    M_after_GA = output_GA.m(end); %  - sim.M_pods;
    [output_1] = NL_interpolator( r_GA , r1 , v_dep_GA , v1 , N_rev1 , TOF1 , M_after_GA ,sim.PS.Isp , sim );
    % 3rd leg - Ast1 -> Ast2
    M_start_2nd_leg = output_1.m(end); %  - sim.M_pods;
    [output_2] = NL_interpolator( r1 , r2 , v1 , v2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );
    % 4th leg - Ast2 -> Ast3
    M_start_3rd_leg = output_2.m(end); %  - sim.M_pods;
    [output_3] = NL_interpolator( r2 , r3 , v2 , v3 , N_rev3 , TOF3 , M_start_3rd_leg ,sim.PS.Isp , sim );
    % 5th leg - Ast3 -> Ast4
    M_start_4th_leg = output_3.m(end); %  - sim.M_pods;
    [output_4] = NL_interpolator( r3 , r4 , v3 , v4 , N_rev4 , TOF4 , M_start_4th_leg ,sim.PS.Isp , sim );
    mass_fract = (output_GA.m(1) - output_4.m(end))/output_GA.m(1);
    % put one after the other, all the thrust profiles
%     T_append = [output_GA.Thrust;
%     output_1.Thrust;
%     output_2.Thrust;
%     output_3.Thrust;
%     output_4.Thrust]; 
%     T = sqrt(T_append(:,1).^2 + T_append(:,3).^2);
% 
%     if abs(max()) <= sim.max_Available_Thrust % bepi colombo is 250 mN
        if mass_fract > 0 && mass_fract < 1 %17 kg of payload
            obj_fun = mass_fract;
            disp('success')
        else
            obj_fun = obj_fun + 10;
        end
%     else
%         obj_fun = obj_fun + 21;
%     end

end

