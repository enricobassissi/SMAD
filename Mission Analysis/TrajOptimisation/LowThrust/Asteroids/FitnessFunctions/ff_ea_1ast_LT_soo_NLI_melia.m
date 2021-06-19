function obj_fun = ff_ea_1ast_LT_soo_NLI_melia(x,sim,data,power_propulsion_data)
% setting the input times
MJD01 = x(1);
TOF1 = x(2);
MJDP1 = MJD01 + TOF1;

% N REV
N_rev = x(3);

% C3 launcher
v_inf_magn = x(4);
az = x(5);
el = x(6);

% chosing which asteroid to visit
IDP = x(7); %index of permutation, the column of the Permutation Matrix of the asteroids
% asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_1 = data.asteroid_names(IDP);
% asteroid_2 = data.PermutationMatrix(IDP,2);
% asteroid_3 = data.PermutationMatrix(IDP,3);
% asteroid_4 = data.PermutationMatrix(IDP,4);

% Computing position and velocity of the planets in that days
% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;

% passage at 1st ast
MJDP1_dim = MJDP1*sim.TU/(3600*24);
[kep_ast_1] = uNEO2(MJDP1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1, v1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
r1 = r1/sim.DU;
v1 = v1/sim.DU*sim.TU;

v_launcher = v_inf_magn*[cos(el)*cos(az); cos(el)*sin(az); sin(el)] ;
v_dep = v_EA + v_launcher;  %since parabolic escape (vinf = 0)
varr = v1; %since we want to rendezvous

%% NLI
tol_TOF = 0.5; penalty_MF_unfeasible = 0;
penalty_T = 0; penalty_TOF = 0;

[output] = NL_interpolator_of( r_EA , r1 , v_dep , varr , N_rev , TOF1 , sim.M ,sim.PS.Isp , sim );

% if max(abs(output.T_magn)) > sim.max_Available_Thrust
%     penalty_T = max(abs(output.T_magn)) - sim.max_Available_Thrust;
% end
if abs(output.t(end) - TOF1) > tol_TOF
    penalty_TOF = abs(output.t(end) - TOF1);
end

mass_fract = (output.m(1) - output.m(end))/output.m(1);

if mass_fract <= 0 || mass_fract >= 1
    penalty_MF_unfeasible = 100;
end

%% thrust available
if ~isnan(output.Href)
    % t in years
    t = output.t*sim.TU/(3600*24*365);
    % trajectory heliocentric
    r3  = [output.r.*cos(output.theta) output.r.*sin(output.theta) output.z];
    R3 = rotate_local2ecplitic(r_EA,r3,sim.n_sol,output.Href);
    % thrust in heliocentric
    Tlocal_transf_orbit  = [-output.Thrust(:,1).*sin(output.theta), ...
        output.Thrust(:,1).*cos(output.theta), output.Thrust(:,3)];
    Thrust_Helio = rotate_local2ecplitic(r_EA,Tlocal_transf_orbit,sim.n_sol,output.Href);
    T_magn_Helio = vecnorm(Thrust_Helio,2,2).*1e3; % mN
    %thrust available
    T_available = available_thrust(t, R3, Thrust_Helio, power_propulsion_data); % mN
    % new penalty on the thrust, if the max(T available -T required > 0 -> penalty
    % so that we give in case higher thrust at beginning or far from the sun because if we 
    % consider all the time as max thrust the one of worst case scenario we limit the 
    % possibilities like if we have a peak at beginning you don't discard the whole mission
    if max(T_magn_Helio - T_available) > 0
        penalty_T = max(T_magn_Helio - T_available);
    end
else
    penalty_T = 10;
end

obj_fun = mass_fract + penalty_MF_unfeasible + 10*penalty_T + penalty_TOF;

end

