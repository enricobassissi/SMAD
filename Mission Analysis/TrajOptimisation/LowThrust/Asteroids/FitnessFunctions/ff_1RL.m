function obj_fun = ff_1RL(x,sim,data)
%{ 
input variable vector
x = [(1) MJD0,
     (2) TOF1,
     (3) TOF2,
     (4) NREV,
     (5) NREV2,
     (6) IDP,
     (7) v_inf_magn,
     (8) alpha (azimuth),
     (9) beta (elevation), 
     (10) CT1, 
     (11) CT2, 
     (12) TOF3,
     (13) NREV3,
     (14) CT3,
     (15) TOF4,
     (16) NREV4]
%}

%% setting the input times
MJD01 = x(1); % departure time from earth
TOF1 = x(2);
MJDA1 = MJD01 + TOF1; % arrival time on ast 1
CT1 = x(10); % coasting time 1
MJDD1 = MJDA1 + CT1; % departure from ast 1
TOF2 = x(3);
MJDA2 = MJDD1 + TOF2; % arrival ast 2
CT2 = x(11); % coasting time 2
MJDD2 = MJDA2 + CT2; % departure on ast 2
TOF3 = x(12);
MJDA3 = MJDD2 + TOF3; % arrival ast 3
CT3 = x(14); % coasting time 3
MJDD3 = MJDA3 + CT3; % departure on ast 3
TOF4 = x(15);
MJDA4 = MJDD3 + TOF4; % arrival ast 4

obj_fun = 10;

max_duration = 12*365*(3600*24)/sim.TU;
if (TOF1+CT1+TOF2+CT2+TOF3+CT3+TOF4) < max_duration % 12 years max mission time 
    
    % N REV1
    N_rev1 = x(4);
    % N REV2
    N_rev2 = x(5);
    % N REV3
    N_rev3 = x(13);
    % N REV4
    N_rev4 = x(16);

    % C3 launcher
    v_inf_magn = x(7);
    az = x(8);
    elev = x(9);

    % chosing which asteroid to visit
    IDP = x(6); %index of permutation, the column of the Permutation Matrix of the asteroids
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

    % arrival at 1st ast
    MJDA1_dim = MJDA1*sim.TU/(3600*24);
    [kep_ast_A1] = uNEO2(MJDA1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
    [rA1, vA1] = sv_from_coe(kep_ast_A1,ksun); % km, km/s
    rA1 = rA1/sim.DU;
    vA1 = vA1/sim.DU*sim.TU;

    % departure from 1st ast
    MJDD1_dim = MJDD1*sim.TU/(3600*24);
    [kep_ast_D1] = uNEO2(MJDD1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
    [rD1, vD1] = sv_from_coe(kep_ast_D1,ksun); % km, km/s
    rD1 = rD1/sim.DU;
    vD1 = vD1/sim.DU*sim.TU;

    % arrival at 2nd ast
    MJDA2_dim = MJDA2*sim.TU/(3600*24);
    [kep_ast_A2] = uNEO2(MJDA2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
    [rA2, vA2] = sv_from_coe(kep_ast_A2,ksun); % km, km/s
    rA2 = rA2/sim.DU;
    vA2 = vA2/sim.DU*sim.TU;

    % departure from 2nd ast
    MJDD2_dim = MJDD2*sim.TU/(3600*24);
    [kep_ast_D2] = uNEO2(MJDD2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
    [rD2, vD2] = sv_from_coe(kep_ast_D2,ksun); % km, km/s
    rD2 = rD2/sim.DU;
    vD2 = vD2/sim.DU*sim.TU;

    % arrival at 3rd ast
    MJDA3_dim = MJDA3*sim.TU/(3600*24);
    [kep_ast_A3] = uNEO2(MJDA3_dim,asteroid_3,data); % [km,-,rad,rad,rad,wrapped rad]
    [rA3, vA3] = sv_from_coe(kep_ast_A3,ksun); % km, km/s
    rA3 = rA3/sim.DU;
    vA3 = vA3/sim.DU*sim.TU;
    
    % departure from 3rd ast
    MJDD3_dim = MJDD3*sim.TU/(3600*24);
    [kep_ast_D3] = uNEO2(MJDD3_dim,asteroid_3,data); % [km,-,rad,rad,rad,wrapped rad]
    [rD3, vD3] = sv_from_coe(kep_ast_D3,ksun); % km, km/s
    rD3 = rD3/sim.DU;
    vD3 = vD3/sim.DU*sim.TU;

    % arrival at 4th ast
    MJDA4_dim = MJDA4*sim.TU/(3600*24);
    [kep_ast_A4] = uNEO2(MJDA4_dim,asteroid_4,data); % [km,-,rad,rad,rad,wrapped rad]
    [rA4, vA4] = sv_from_coe(kep_ast_A4,ksun); % km, km/s
    rA4 = rA4/sim.DU;
    vA4 = vA4/sim.DU*sim.TU;

    %% Launcher departure variable
%     z = r .* sin(elev);
%     rcoselev = r .* cos(elev);
%     x = rcoselev .* cos(az);
%     y = rcoselev .* sin(az);
    v_launcher = v_inf_magn*[cos(elev)*cos(az); cos(elev)*sin(az); sin(elev)];
    v_dep = v_EA + v_launcher;  %if parabolic escape (v_extra = 0)
   
    %% NLI
    % 1st leg - Earth -> Ast 1
    [output_1] = NL_interpolator( r_EA , rA1 , v_dep , vA1 , N_rev1 , TOF1 , sim.M ,sim.PS.Isp , sim );
    sim.M_pods = 0;
    if ~isnan(output_1.Thrust) % if is not nan -> it's a valid solution
        % 2nd leg - Ast1 -> Ast2
        M_start_2nd_leg = output_1.m(end) - sim.M_pods;
        [output_2] = NL_interpolator( rD1 , rA2 , vD1 , vA2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );

        if ~isnan(output_2.Thrust) % if is not nan -> it's a valid solution
            % 3rd leg - Ast2 -> Ast3
            M_start_3rd_leg = output_2.m(end) - sim.M_pods;
            [output_3] = NL_interpolator( rD2 , rA3 , vD2 , vA3 , N_rev3 , TOF3 , M_start_3rd_leg ,sim.PS.Isp , sim );

            if ~isnan(output_3.Thrust) % if is not nan -> it's a valid solution
                % 4th leg - Ast3 -> Ast4
                M_start_4th_leg = output_3.m(end) - sim.M_pods;
                [output_4] = NL_interpolator( rD3 , rA4 , vD3 , vA4 , N_rev4 , TOF4 , M_start_4th_leg ,sim.PS.Isp , sim );
                
                if ~isnan(output_4.Thrust) % if is not nan -> it's a valid solution
                    mass_fract = (output_1.m(1) - output_4.m(end))/output_1.m(1);

                    % put one after the other, all the thrust profiles
                    T_append = [output_1.Thrust(:,1),output_1.Thrust(:,2),output_1.Thrust(:,3);
                                output_2.Thrust(:,1),output_2.Thrust(:,2),output_2.Thrust(:,3);
                                output_3.Thrust(:,1),output_3.Thrust(:,2),output_3.Thrust(:,3);
                                output_4.Thrust(:,1),output_4.Thrust(:,2),output_4.Thrust(:,3)]; 
                    T = sqrt(T_append(:,1).^2 + T_append(:,3).^2);

                    if abs(max(T)) <= 0.1 % bepi colombo is 250 mN
                        if mass_fract > 0 && mass_fract < 1 %17 kg of payload
                            obj_fun = mass_fract;
                            disp('success')
                        else
%                             obj_fun = 1e4; % exception of mass fraction wrong
%                             obj_fun = 10+1e3*(mass_fract-0.6)^2; % exception of mass fraction wrong
%                             disp('mass_fract error')
                            obj_fun = obj_fun + 10;
                        end
                    else
%                         obj_fun = 1e4; % exception of max thrust exceeded
%                         obj_fun = 10+1e3*(abs(max(T))-0.015)^2;
%                         disp('max T error')
                        obj_fun = obj_fun + 21;
                    end
                else
%                     disp('4th leg error') % exception of 4th leg wrong
                    obj_fun = obj_fun + 32;
                end
            else
%                 disp('3rd leg error') % exception of 3rd leg wrong
                obj_fun = obj_fun + 43;
            end
        else
%             disp('2nd leg error') % exception of 2nd leg wrong
            obj_fun = obj_fun + 54;
        end
    else
%         disp('1st leg error') % exception of 1st leg wrong
        obj_fun = obj_fun + 65;
    end
else
%     obj_fun = 1e4; % exception of too long mission
%     obj_fun = 1e2*(TOF1+CT1+TOF2+CT2+TOF3 - max_duration)^2;
    disp('max mission duration error')
    obj_fun = obj_fun + 76;
end

end

