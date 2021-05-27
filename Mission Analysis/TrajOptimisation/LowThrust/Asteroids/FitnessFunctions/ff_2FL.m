function obj_fun = ff_2FL(x,sim,data)

%% setting the input times
MJD01 = x(1); % departure time from earth
TOF1 = x(2);
MJDP1 = MJD01 + TOF1; % passage time on ast 1
TOF2 = x(3);
MJDP2 = MJDP1 + TOF2; % passage ast 2

TOFa = x(4);
MJDPa = MJD01 + TOFa; % passage ast 3
TOFb = x(5);
MJDPb = MJDPa + TOFb; % passage ast 4

obj_fun = 100; %%

max_duration = 12*365*(3600*24)/sim.TU;
penalty_MAX_DURATION = 0;
if max(TOF1+TOF2,TOFa+TOFb) > max_duration
    penalty_MAX_DURATION = max(TOF1+TOF2,TOFa+TOFb) - max_duration; % 12 years max mission time 
end
    
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
el = x(14);


% chosing which asteroids to visit

% 1st SPACECRAFT ASTEROIDS OBJECTIVES
IDP = x(10); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_2 = data.PermutationMatrix(IDP,2);

% 2nd SPACECRAFT ASTEROIDS OBJECTIVES
IDP_temp_2 = x(11); % index for 2nd permutation matrix to be built inside depending on the first 2 selected asteroids
asteroid_sequence = [asteroid_1,asteroid_2];
TF = contains(data.asteroid_names,asteroid_sequence);
data_elements_matrix_2SC = data.data_elements_matrix(~TF,:);
[~, PermutationMatrix_2SC, HowMany_2SC] = ...
            sequences_local_pruning2(data_elements_matrix_2SC, data.p_number);
IDP2 = ceil(IDP_temp_2*HowMany_2SC/100);
asteroid_a = PermutationMatrix_2SC(IDP2,1);
asteroid_b = PermutationMatrix_2SC(IDP2,2);


% asteroid1 flyby
v_inf_ast1_magn = x(15);
az_ast1 = x(16);
el_ast1 = x(17);

% asteroid2 flyby
v_inf_ast2_magn = x(18);
az_ast2 = x(19);
el_ast2 = x(20);

% asteroida flyby
v_inf_asta_magn = x(21);
az_asta = x(22);
el_asta = x(23);

% asteroidb flyby
v_inf_astb_magn = x(24);
az_astb = x(25);
el_astb = x(26);
     
    
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

% passage at a ast
MJDPa_dim = MJDPa*sim.TU/(3600*24);
[kep_ast_a] = uNEO2(MJDPa_dim,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[ra, va] = sv_from_coe(kep_ast_a,ksun); % km, km/s
ra = ra/sim.DU;
va = va/sim.DU*sim.TU;

% passage at b ast
MJDPb_dim = MJDPb*sim.TU/(3600*24);
[kep_ast_b] = uNEO2(MJDPb_dim,asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
[rb, vb] = sv_from_coe(kep_ast_b,ksun); % km, km/s
rb = rb/sim.DU;
vb = vb/sim.DU*sim.TU;


%% Relative velocity -> heliocentric velocity
% Launcher departure variable
v_extralaunch = v_inf_magn*[cos(el)*cos(az); cos(el)*sin(az); sin(el)];
v_dep = v_EA + v_extralaunch;  %if parabolic escape (v_extralaunch = 0)

% Flyby ast 1
v_rel_ast1 = v_inf_ast1_magn* [cos(el_ast1)*cos(az_ast1); cos(el_ast1)*sin(az_ast1); sin(el_ast1)];
v_abs_ast1 = v1 + v_rel_ast1;

% Flyby ast 2
v_rel_ast2 = v_inf_ast2_magn* [cos(el_ast2)*cos(az_ast2); cos(el_ast2)*sin(az_ast2); sin(el_ast2)];
v_abs_ast2 = v2 + v_rel_ast2;

% Flyby ast a
v_rel_asta = v_inf_asta_magn* [cos(el_asta)*cos(az_asta); cos(el_asta)*sin(az_asta); sin(el_asta)];
v_abs_asta = va + v_rel_asta;

% Flyby ast b
v_rel_astb = v_inf_astb_magn* [cos(el_astb)*cos(az_astb); cos(el_astb)*sin(az_astb); sin(el_astb)];
v_abs_astb = vb + v_rel_astb;

%% NLI
tol_TOF = 1; % 1 TU means approx 60 days
penalty_T_leg1 = 0; penalty_T_leg2 = 0; penalty_T_lega = 0; penalty_T_legb = 0; 
penalty_TOF_leg1 = 0; penalty_TOF_leg2 = 0; penalty_TOF_lega = 0; penalty_TOF_legb = 0;
% SC1
% 1st leg - Earth -> ast1
[output_1] = NL_interpolator_of( r_EA , r1 , v_dep , v_abs_ast1, N_rev1 , TOF1 , sim.M1 ,sim.PS.Isp , sim );
if max(abs(output_1.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg1 = abs(max(output_1.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_1.t(end) - TOF1) > tol_TOF
    penalty_TOF_leg1 = abs(output_1.t(end) - TOF1);
end
% if ~isnan(output_1.Thrust)
% 2nd leg - Ast1 -> Ast2
M_start_2nd_leg = output_1.m(end); %  - sim.M_pods;
[output_2] = NL_interpolator_of( r1 , r2 , v_abs_ast1 , v_abs_ast2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );
if max(abs(output_2.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg2 = abs(max(output_2.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_2.t(end) - TOF2) > tol_TOF
    penalty_TOF_leg2 = abs(output_2.t(end) - TOF2);
end

% SC2
% 1st leg - Earth -> asta
[output_a] = NL_interpolator_of( r_EA , ra , v_dep , v_abs_asta, N_reva , TOFa , sim.M2 ,sim.PS.Isp , sim );
if max(abs(output_a.T_magn)) > sim.max_Available_Thrust
    penalty_T_lega = 10*abs(max(output_a.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_a.t(end) - TOFa) > tol_TOF
    penalty_TOF_lega = abs(output_a.t(end) - TOFa);
end
% if ~isnan(output_a.Thrust)
% 2nd leg - Asta -> Astb
M_start_b_leg = output_a.m(end); %  - sim.M_pods;
[output_b] = NL_interpolator_of( ra , rb , v_abs_asta , v_abs_astb , N_revb , TOFb , M_start_b_leg ,sim.PS.Isp , sim );
if max(abs(output_b.T_magn)) > sim.max_Available_Thrust
    penalty_T_legb = 10*abs(max(output_b.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_b.t(end) - TOFb) > tol_TOF
    penalty_TOF_legb = abs(output_b.t(end) - TOFb);
end


mass_fract_SC1 = (output_1.m(1) - output_2.m(end))/output_1.m(1);
mass_fract_SC2 = (output_a.m(1) - output_b.m(end))/output_a.m(1);

penalty_MF_feasibility = 0;
if mass_fract_SC1 < 0 || mass_fract_SC1 > 1 || mass_fract_SC2 < 0 || mass_fract_SC2 > 1
    penalty_MF_feasibility = 100;
end
%         if mass_fract_SC2 > 0 && mass_fract_SC2 < 1 
avg_mass_fraction = (mass_fract_SC1+mass_fract_SC2)/2;
MF = max(mass_fract_SC1,mass_fract_SC2) + abs(mass_fract_SC1 - avg_mass_fraction) + ...
        abs(mass_fract_SC2 - avg_mass_fraction); % cosi sono piu o meno uguali


obj_fun = MF + penalty_MAX_DURATION + penalty_MF_feasibility + ...
    10*(penalty_T_leg1 + penalty_T_leg2 + penalty_T_lega + penalty_T_legb) + ...
    penalty_TOF_leg1 + penalty_TOF_leg2 + penalty_TOF_lega + penalty_TOF_legb;

    
end



                 
                        
