function obj_fun = ff_1FL_GA_Ale(x,sim,data)

%% setting the input times
MJD01 = x(1); % departure time from earth
TOFGA = x(14);
MJDPGA = MJD01 + TOFGA; % passage time GA
TOF1 = x(2);
MJDP1 = MJDPGA + TOF1; % passage time on ast 1
TOF2 = x(3);
MJDP2 = MJDP1 + TOF2; % passage ast 2
TOF3 = x(4);
MJDP3 = MJDP2 + TOF3; % passage ast 3
TOF4 = x(5);
MJDP4 = MJDP3 + TOF4; % passage ast 4

obj_fun = 1000; %%

%max_duration = 12*365*(3600*24)/sim.TU;
%if (TOFGA+TOF1+TOF2+TOF3+TOF4) < max_duration % 12 years max mission time 
    
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
el = x(13);

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

% asteroid1 flyby
v_inf_ast1_magn = x(21);
az_ast1 = x(22);
el_ast1 = x(23);

% asteroid2 flyby
v_inf_ast2_magn = x(24);
az_ast2 = x(25);
el_ast2 = x(26);

% asteroid3 flyby
v_inf_ast3_magn = x(27);
az_ast3 = x(28);
el_ast3 = x(29);

% asteroid4 flyby
v_inf_ast4_magn = x(30);
az_ast4 = x(31);
el_ast4 = x(32);
     
    
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


%% Relative velocity -> heliocentric velocity
% Launcher departure variable
v_extralaunch = v_inf_magn*[cos(el)*cos(az); cos(el)*sin(az); sin(el)];
v_dep = v_EA + v_extralaunch;  %if parabolic escape (v_extralaunch = 0)

% Gravity assist
v_arr_GA = v_GA + v_inf_magn2*[cos(el2)*cos(az2); cos(el2)*sin(az2); sin(el2)];
v_dep_GA = v_GA + v_inf_magn2*[cos(el3)*cos(az3); cos(el3)*sin(az3); sin(el3)];

% Flyby ast 1
v_rel_ast1 = v_inf_ast1_magn* [cos(el_ast1)*cos(az_ast1); cos(el_ast1)*sin(az_ast1); sin(el_ast1)];
v_abs_ast1 = v1 + v_rel_ast1;

% Flyby ast 2
v_rel_ast2 = v_inf_ast2_magn* [cos(el_ast2)*cos(az_ast2); cos(el_ast2)*sin(az_ast2); sin(el_ast2)];
v_abs_ast2 = v2 + v_rel_ast2;

% Flyby ast 3
v_rel_ast3 = v_inf_ast3_magn* [cos(el_ast3)*cos(az_ast3); cos(el_ast3)*sin(az_ast3); sin(el_ast3)];
v_abs_ast3 = v3 + v_rel_ast3;

% Flyby ast 4
v_rel_ast4 = v_inf_ast4_magn* [cos(el_ast4)*cos(az_ast4); cos(el_ast4)*sin(az_ast4); sin(el_ast4)];
v_abs_ast4 = v4 + v_rel_ast4;

%% NLI
% 1st leg - Earth -> GA
[output_GA] = NL_interpolator( r_EA , r_GA , v_dep , v_arr_GA , N_revGA , TOFGA , sim.M ,sim.PS.Isp , sim );
    
if ~isnan(output_GA.Thrust) 
        % 2nd leg - GA -> Ast 1
        if sim.ID_FLYBY == 3
            RPlanet_flyby = astroConstants(23); % Radius_Earth, km
            muPlanet_flyby = astroConstants(13); % muEarth, km^3/s^2
            R_lim_from_planet = 500; % km, for earth is ok to avoid atmosphere
            R_SOI_PL = 0.929*1e6; % km
        elseif sim.ID_FLYBY == 4
            RPlanet_flyby = astroConstants(24); % Radius_mars, km
            muPlanet_flyby = astroConstants(14); % mu mars, km^3/s^2
            R_lim_from_planet = 200; % km, for mars is ok to avoid atmosphere
            R_SOI_PL = 0.578*1e6; %km
        end
        v_arr_GA_dim = v_arr_GA.*sim.DU./sim.TU;
        v_dep_GA_dim = v_dep_GA.*sim.DU./sim.TU;
        [delta_v_p,rp] = flyby(RPlanet_flyby, muPlanet_flyby,R_lim_from_planet, ...
                          MJDPGA_dim, v_arr_GA_dim, v_dep_GA_dim, sim.ID_FLYBY,R_SOI_PL);

                      
        if ~strcmp(string(delta_v_p), 'Not found')
            M_after_GA = output_GA.m(end);
            [output_1] = NL_interpolator( r_GA , r1 , v_dep_GA , v_abs_ast1, N_rev1 , TOF1 , M_after_GA ,sim.PS.Isp , sim );
       
                if ~isnan(output_1.Thrust)
                        % 3rd leg - Ast1 -> Ast2
                        M_start_2nd_leg = output_1.m(end); %  - sim.M_pods;
                        [output_2] = NL_interpolator( r1 , r2 , v_abs_ast1 , v_abs_ast2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );

                        if ~isnan(output_2.Thrust) 
                            % 4th leg - Ast2 -> Ast3
                            M_start_3rd_leg = output_2.m(end); %  - sim.M_pods;
                            [output_3] = NL_interpolator( r2 , r3 , v_abs_ast2 , v_abs_ast3 , N_rev3 , TOF3 , M_start_3rd_leg ,sim.PS.Isp , sim );

                            if ~isnan(output_3.Thrust) 
                                % 5th leg - Ast3 -> Ast4
                                M_start_4th_leg = output_3.m(end); %  - sim.M_pods;
                                [output_4] = NL_interpolator( r3 , r4 , v_abs_ast3 , v_abs_ast4 , N_rev4 , TOF4 , M_start_4th_leg ,sim.PS.Isp , sim );

                                if ~isnan(output_4.Thrust) 
                                   mass_fract = (output_GA.m(1) - output_4.m(end))/output_GA.m(1);
                                    

                                    % put one after the other, all the thrust profiles
                                    T_append = [output_GA.Thrust(:,1),output_GA.Thrust(:,2),output_GA.Thrust(:,3);
                                                output_1.Thrust(:,1),output_1.Thrust(:,2),output_1.Thrust(:,3);
                                                output_2.Thrust(:,1),output_2.Thrust(:,2),output_2.Thrust(:,3);
                                                output_3.Thrust(:,1),output_3.Thrust(:,2),output_3.Thrust(:,3);
                                                output_4.Thrust(:,1),output_4.Thrust(:,2),output_4.Thrust(:,3)];

                                            
                                     T = sqrt(T_append(:,1).^2 + T_append(:,3).^2);

                                        if mass_fract > 0 && mass_fract < 1 %17 kg of payload
                                            %obj_fun = mass_fract;
                                            obj_fun = mass_fract + sim.c*(abs(max(T))-sim.maxT)^2;
                                            %disp('success')
                                        else
                                            obj_fun = obj_fun + 11;
                                        end
                                else
                                    obj_fun = obj_fun + 32;
                                end
                            else
                                obj_fun = obj_fun + 43;
                            end
                        else
                            obj_fun = obj_fun + 54;
                        end
                else
                    obj_fun = obj_fun + 65;
                end
        else
            obj_fun = obj_fun + 76;
        end
else
        obj_fun = obj_fun + 87;
end
        
    

end
   
