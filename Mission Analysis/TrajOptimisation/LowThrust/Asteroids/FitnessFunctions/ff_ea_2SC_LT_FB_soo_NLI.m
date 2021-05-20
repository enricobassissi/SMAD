function obj_fun = ff_ea_2SC_LT_FB_soo_NLI(x,sim,data)
%{ 
input variable vector
x = [(1) MJD0,
     (2) TOF1,
     (3) TOF2,
     (4) TOFa,
     (5) TOFb,
     (6) NREV,
     (7) NREV2,
     (8) NREVa,
     (9) NREVb,
     (10) IDP,
     (11) IDP2_0_100,
     (12) v_inf_magn,
     (13) azimuth,
     (14) elevation]
%}

% Nomenclature
% quantities for SpaceCraft1 will have numbers (1,2,...)
% quantities for SC 2 will have letters (a,b,...)

% setting the input times
MJD01 = x(1); % departure time for both the sc

% 1st spacecraft characteristic times
TOF1 = x(2); % tof sc1 to 1st asteroid
MJDP1 = MJD01 + TOF1; % mjd2000 passage of 1st sc on ast 1
TOF2 = x(3); % tof sc1 to 2nd asteroid
MJDP2 = MJDP1 + TOF2; % mjd2000 passage of 1st sc on ast 2

% 2nd spacecraft characteristic times
TOFa = x(4); % tof sc2 to 1st asteroid
MJDPa = MJD01 + TOFa; 
TOFb =  x(5); % tof sc2 to 2nd asteroid
MJDPb = MJDPa + TOFb; 

obj_fun = 10;

max_duration = 12*365*(3600*24)/sim.TU;
if max(TOF1+TOF2,TOFa+TOFb) < max_duration % 12 years max mission time 
    
    % N REV1
    N_rev1 = x(6);
    % N REV2
    N_rev2 = x(7);
    % N REVa
    N_reva = x(8);
    % N REVb
    N_revb = x(9);

    % C3 launcher
    v_inf_magn = x(12);
    az = x(13);
    elev = x(14);

    %% choosing which asteroid to visit
    % 1ST SPACECRAFT ASTEROID OBJECTIVES
    IDP1 = x(10); %index of permutation, the column of the Permutation Matrix of the asteroids
    asteroid_1 = data.PermutationMatrix(IDP1,1);
    asteroid_2 = data.PermutationMatrix(IDP1,2);

    % 2ND SPACECRAFT ASTEROID OBJECTIVES
    IDP_temp_2 = x(11); % index for 2nd permutation matrix to be built inside depending on the first 2 selected asteroids
    asteroid_sequence = [asteroid_1,asteroid_2];
    TF = contains(data.asteroid_names,asteroid_sequence);
    data_elements_matrix_2SC = data.data_element_matrix(~TF,:);
    [~, PermutationMatrix_2SC, HowMany_2SC] = ...
                sequences_local_pruning(data_elements_matrix_2SC, data.p_number);
    IDP2 = ceil(IDP_temp_2*HowMany_2SC/100);
    asteroid_a = PermutationMatrix_2SC(IDP2,1);
    asteroid_b = PermutationMatrix_2SC(IDP2,2);

    %% Computing position and velocity of the planets in that days
    % Departure from Earth
    MJD01_dim = MJD01*sim.TU/(3600*24);
    [kep_EA,ksun] = uplanet(MJD01_dim, 3);
    [r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
    r_EA = r_EA/sim.DU;
    v_EA = v_EA/sim.DU*sim.TU;

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

    % passage at a_th ast
    MJDPa_dim = MJDPa*sim.TU/(3600*24);
    [kep_ast_a] = uNEO2(MJDPa_dim,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
    [ra, va] = sv_from_coe(kep_ast_a,ksun); % km, km/s
    ra = ra/sim.DU;
    va = va/sim.DU*sim.TU;

    % passage at b_th ast
    MJDPb_dim = MJDPb*sim.TU/(3600*24);
    [kep_ast_b] = uNEO2(MJDPb_dim,asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
    [rb, vb] = sv_from_coe(kep_ast_b,ksun); % km, km/s
    rb = rb/sim.DU;
    vb = vb/sim.DU*sim.TU;

    %% Launcher departure variable
%     z = r .* sin(elev);
%     rcoselev = r .* cos(elev);
%     x = rcoselev .* cos(az);
%     y = rcoselev .* sin(az);
    v_launcher = v_inf_magn*[cos(elev)*cos(az); cos(elev)*sin(az); sin(elev)];
    v_dep = v_EA + v_launcher;  %if parabolic escape (v_extra = 0)
   
    %% NLI
    % SC1 
    % 1st leg - Earth -> Ast 1
    [output_1] = NL_interpolator( r_EA , r1 , v_dep , v1 , N_rev1 , TOF1 , sim.M1 ,sim.PS.Isp , sim );
    if ~isnan(output_1.Thrust) % if is not nan -> it's a valid solution
        % 2nd leg - Ast1 -> Ast2
        M_start_2nd_leg = output_1.m(end); %  - sim.M_pods
        [output_2] = NL_interpolator( r1 , r2 , v1 , v2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );
        if ~isnan(output_2.Thrust) % if is not nan -> it's a valid solution
            % SC 2
            [output_a] = NL_interpolator( r_EA , ra , v_dep , va , N_reva , TOFa , sim.M2 ,sim.PS.Isp , sim );
            if ~isnan(output_a.Thrust) % if is not nan -> it's a valid solution
                % a_th leg - Earth -> Ast_a
                M_start_b_th_leg = output_a.m(end); %  - sim.M_pods
                [output_b] = NL_interpolator( ra , rb , va , vb , N_revb , TOFb , M_start_b_th_leg ,sim.PS.Isp , sim );
                if ~isnan(output_b.Thrust) % if is not nan -> it's a valid solution
                    % mass fractions
                    mass_fract_SC1 = (output_1.m(1) - output_2.m(end))/output_1.m(1);
                    mass_fract_SC2 = (output_a.m(1) - output_b.m(end))/output_a.m(1);
                    % thrust profiles
                    T_append_SC1 = [output_1.Thrust(:,1),output_1.Thrust(:,2),output_1.Thrust(:,3);
                                    output_2.Thrust(:,1),output_2.Thrust(:,2),output_2.Thrust(:,3)]; 
                    T_SC1 = sqrt(T_append_SC1(:,1).^2 + T_append_SC1(:,3).^2);
                    
                    T_append_SC2 = [output_a.Thrust(:,1),output_a.Thrust(:,2),output_a.Thrust(:,3);
                                    output_b.Thrust(:,1),output_b.Thrust(:,2),output_b.Thrust(:,3)]; 
                    T_SC2 = sqrt(T_append_SC2(:,1).^2 + T_append_SC2(:,3).^2);

                    if abs(max(T_SC1)) <= 0.05 % bepi colombo is 250 mN
                        if abs(max(T_SC2)) <= 0.05
                            if mass_fract_SC1 > 0 && mass_fract_SC1 < 1 %17 kg of payload
                                if mass_fract_SC2 > 0 && mass_fract_SC2 < 1 %17 kg of payload
                                    obj_fun = max(mass_fract_SC1,mass_fract_SC2); % + mean(mf_sc1,mf_sc2), cosi sono piu o meno uguali
                                    disp('success')
                                else
                                    obj_fun = obj_fun + 10; % error in the mass fraction of SC2
                                end
                            else
                                obj_fun = obj_fun + 11; % error in the mass fraction of SC1
                            end
                        else
                            obj_fun = obj_fun + 22; % error in the max thrust of SC2
                        end
                    else
                        obj_fun = obj_fun + 23; % error in the max thrust of SC1
                    end
                else
                    obj_fun = obj_fun + 34; % error in the 2nd leg of SC2
                end
            else
                obj_fun = obj_fun + 35; % error in the 1st leg of SC2
            end
        else
            obj_fun = obj_fun + 46; % error in the 2nd leg of SC1
        end
    else
        obj_fun = obj_fun + 47; % error in the 1st leg of SC1
    end
else
    disp('max mission duration error')
    obj_fun = obj_fun + 58;
end

end

