function [sol] = plot_ff_2RL_all_indietro_AC2_flex(x,sim,data)
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
     (10) v_inf_magn,
     (11) azimuth,
     (12) elevation,
     (13) CT1,
     (14) CT2]
%}

% Nomenclature
% quantities for SpaceCraft1 will have numbers (1,2,...)
% quantities for SC 2 will have letters (a,b,...)

% setting the input times
MJD01 = x(1); % departure time for both the sc

% 1st spacecraft characteristic times
TOF1 = x(2); % tof sc1 to 1st asteroid
MJDA1 = MJD01 + TOF1; % mjd2000 passage of 1st sc on ast 1
CT1 = x(13);
MJDD1 = MJDA1 + CT1;
TOF2 = x(3); % tof sc1 to 2nd asteroid
MJDA2 = MJDD1 + TOF2; % mjd2000 passage of 1st sc on ast 2

% 2nd spacecraft characteristic times
TOFa = x(4); % tof sc2 to 1st asteroid
MJDAa = MJD01 + TOFa; 
CTa = x(14);
MJDDa = MJDAa + CTa;
TOFb =  x(5); % tof sc2 to 2nd asteroid
MJDAb = MJDDa + TOFb; 

max_duration = 12*365*(3600*24)/sim.TU;
penalty_MAX_DURATION = 0;
if max(TOF1+CT1+TOF2,TOFa+CTa+TOFb) > max_duration
    penalty_MAX_DURATION = max(TOF1+CT1+TOF2,TOFa+CTa+TOFb) - max_duration; % 12 years max mission time 
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
v_inf_magn = x(10);
az = x(11);
elev = x(12);

%% choosing which asteroid to visit
if data.IDX_WhichSpacecraftWeAreChanging == 1
    if data.IDX_WhichAstWeAreChanging == 1
        % 1st sc ast change
        asteroid_1 = data.AstToKeep{1,1};
        asteroid_2 = data.NextAstSameSpacecraft;
    elseif data.IDX_WhichAstWeAreChanging == 2
        % 1st sc ast change
        asteroid_1 = data.NextAstSameSpacecraft;
        asteroid_2 = data.AstToKeep{1,2};
    end
    % 2nd sc's asteroid are defined
    asteroid_a = data.OtherSpacecraftFixedAst(1);
    asteroid_b = data.OtherSpacecraftFixedAst(2);
    
elseif data.IDX_WhichSpacecraftWeAreChanging == 2
    if data.IDX_WhichAstWeAreChanging == 1
        % 2nd sc's asteroid are changing
        asteroid_a = data.AstToKeep{2,1};
        asteroid_b = data.NextAstSameSpacecraft;
    elseif data.IDX_WhichAstWeAreChanging == 2
        % 2nd sc's asteroid are changing
        asteroid_a = data.NextAstSameSpacecraft;
        asteroid_b = data.AstToKeep{2,2};
    end
    % 1st sc ast fixed
    asteroid_1 = data.OtherSpacecraftFixedAst(1);
    asteroid_2 = data.OtherSpacecraftFixedAst(2);
end

%% Computing position and velocity of the planets in that days
% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;

% ARRIVAL at 1st ast
MJDA1_dim = MJDA1*sim.TU/(3600*24);
[kep_ast_A1] = uNEO3(MJDA1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[rA1, vA1] = sv_from_coe(kep_ast_A1,ksun); % km, km/s
rA1 = rA1/sim.DU;
vA1 = vA1/sim.DU*sim.TU;
% DEPARTURE at 1st ast
MJDD1_dim = MJDD1*sim.TU/(3600*24);
[kep_ast_D1] = uNEO3(MJDD1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[rD1, vD1] = sv_from_coe(kep_ast_D1,ksun); % km, km/s
rD1 = rD1/sim.DU;
vD1 = vD1/sim.DU*sim.TU;

% ARRIVAL at 2nd ast
MJDA2_dim = MJDA2*sim.TU/(3600*24);
[kep_ast_A2] = uNEO3(MJDA2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[rA2, vA2] = sv_from_coe(kep_ast_A2,ksun); % km, km/s
rA2 = rA2/sim.DU;
vA2 = vA2/sim.DU*sim.TU;

% ARRIVAL at a_th ast
MJDAa_dim = MJDAa*sim.TU/(3600*24);
[kep_ast_Aa] = uNEO3(MJDAa_dim,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[rAa, vAa] = sv_from_coe(kep_ast_Aa,ksun); % km, km/s
rAa = rAa/sim.DU;
vAa = vAa/sim.DU*sim.TU;
% DEPARTURE at a_th ast
MJDDa_dim = MJDDa*sim.TU/(3600*24);
[kep_ast_Da] = uNEO3(MJDDa_dim,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[rDa, vDa] = sv_from_coe(kep_ast_Da,ksun); % km, km/s
rDa = rDa/sim.DU;
vDa = vDa/sim.DU*sim.TU;

% passage at b_th ast
MJDAb_dim = MJDAb*sim.TU/(3600*24);
[kep_ast_Ab] = uNEO3(MJDAb_dim,asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
[rAb, vAb] = sv_from_coe(kep_ast_Ab,ksun); % km, km/s
rAb = rAb/sim.DU;
vAb = vAb/sim.DU*sim.TU;

%% Launcher departure variable
v_launcher = v_inf_magn*[cos(elev)*cos(az); cos(elev)*sin(az); sin(elev)];
v_dep = v_EA + v_launcher;  %if parabolic escape (v_extra = 0)

%% NLI
tol_TOF = 1; % 1 TU means approx 60 days
penalty_T_leg1 = 0; penalty_T_leg2 = 0; penalty_T_lega = 0; penalty_T_legb = 0; 
penalty_TOF_leg1 = 0; penalty_TOF_leg2 = 0; penalty_TOF_lega = 0; penalty_TOF_legb = 0;

% ------------ SC1 ------------- %
%   2nd leg - Ast1 -> Ast2
[output_2] = NL_interpolator_of( rD1 , rA2 , vD1 , vA2 , N_rev2 , TOF2 , sim.M1_end ,sim.PS.Isp , sim );
if max(abs(output_2.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg2 = max(abs(output_2.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_2.t(end) - TOF2) > tol_TOF
    penalty_TOF_leg2 = abs(output_2.t(end) - TOF2);
end

% 1st leg - Earth -> Ast 1
M_start_1st_leg = output_2.m(1) + sim.M_pods;
[output_1] = NL_interpolator_of( r_EA , rA1 , v_dep , vA1 , N_rev1 , TOF1 , M_start_1st_leg ,sim.PS.Isp , sim );
if max(abs(output_1.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg1 = max(abs(output_1.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_1.t(end) - TOF1) > tol_TOF
    penalty_TOF_leg1 = abs(output_1.t(end) - TOF1);
end

% ------ SC 2 ------ %
% bth leg - Ast a -> Ast b
[output_b] = NL_interpolator_of( rDa , rAb , vDa , vAb , N_revb , TOFb , sim.M2_end ,sim.PS.Isp , sim );
if max(abs(output_b.T_magn)) > sim.max_Available_Thrust
    penalty_T_legb = 10*max(abs(output_b.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_b.t(end) - TOFb) > tol_TOF
    penalty_TOF_legb = abs(output_b.t(end) - TOFb);
end

% a_th leg - Earth -> Ast_a
M_start_a_th_leg = output_b.m(1) + sim.M_pods;
[output_a] = NL_interpolator_of( r_EA , rAa , v_dep , vAa , N_reva , TOFa , M_start_a_th_leg ,sim.PS.Isp , sim );
if max(abs(output_a.T_magn)) > sim.max_Available_Thrust
    penalty_T_lega = 10*max(abs(output_a.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_a.t(end) - TOFa) > tol_TOF
    penalty_TOF_lega = abs(output_a.t(end) - TOFa);
end

%% mass fractions
g0 = 9.81;
Isp = sim.PS.Isp*sim.TU;

% --- SC1
sol.mass_depleted_Leg2 = output_2.m(1) - output_2.m(end);
sol.dV_associated_Leg2 = -g0*Isp*log(output_2.m(end)/output_2.m(1)); % -ve*ln(m_final/m_initial)

sol.mass_depleted_Leg1 = output_1.m(1) - output_1.m(end);
sol.dV_associated_Leg1 = -g0*Isp*log(output_1.m(end)/output_1.m(1)); % -ve*ln(m_final/m_initial)

sol.tot_mass_depleted_SC1 = sol.mass_depleted_Leg1+sol.mass_depleted_Leg2;
sol.mass_fract_SC1 = sol.tot_mass_depleted_SC1/output_1.m(1);

% --- SC2
sol.mass_depleted_Lega = output_a.m(1) - output_a.m(end);
sol.dV_associated_Lega = -g0*Isp*log(output_a.m(end)/output_a.m(1)); % -ve*ln(m_final/m_initial)

sol.mass_depleted_Legb = output_b.m(1) - output_b.m(end);
sol.dV_associated_Legb = -g0*Isp*log(output_b.m(end)/output_b.m(1)); % -ve*ln(m_final/m_initial)

sol.tot_mass_depleted_SC2 = sol.mass_depleted_Lega+sol.mass_depleted_Legb;
sol.mass_fract_SC2 = sol.tot_mass_depleted_SC2/output_a.m(1);

%% obj fun
Penalty_about_feasibility_of_mass_fraction = 0;
if sol.mass_fract_SC1 < 0 || sol.mass_fract_SC1 > 1 || sol.mass_fract_SC2 < 0 || sol.mass_fract_SC2 > 1 
    Penalty_about_feasibility_of_mass_fraction = 100;
end

avg_mass_fraction = (sol.mass_fract_SC1+ sol.mass_fract_SC2)/2;
MF = max(sol.mass_fract_SC1,sol.mass_fract_SC2) + 10*abs(sol.mass_fract_SC1 - avg_mass_fraction) + ...
        10*abs(sol.mass_fract_SC2 - avg_mass_fraction); % cosi sono piu o meno uguali

sol.obj_fun(1) = MF + penalty_MAX_DURATION + Penalty_about_feasibility_of_mass_fraction;

sol.obj_fun(2) = 100*(penalty_T_leg1 + penalty_T_leg2 + penalty_T_lega + penalty_T_legb) + ...
    penalty_MAX_DURATION + penalty_TOF_leg1 + penalty_TOF_leg2 + penalty_TOF_lega + penalty_TOF_legb;

%% -- Thrust
sol.T_1 = output_1.Thrust;
sol.T_1_magn  = sqrt(sol.T_1(:,1).^2 + sol.T_1(:,3).^2);
sol.T_2 = output_2.Thrust;
sol.T_2_magn  = sqrt(sol.T_2(:,1).^2 + sol.T_2(:,3).^2);
sol.T_a = output_a.Thrust;
sol.T_a_magn  = sqrt(sol.T_a(:,1).^2 + sol.T_a(:,3).^2);
sol.T_b = output_b.Thrust;
sol.T_b_magn  = sqrt(sol.T_b(:,1).^2 + sol.T_b(:,3).^2);

end

