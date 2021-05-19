function [output,r_encounter,v_encounter] = plot_ff_ea_1ast_LT_GA_FB_soo_NLI(x,sim,data)
% setting the input times
MJD01 = x(1);
TOFGA = x(2);
MJDPGA = MJD01 + TOFGA;
TOF1 = x(3);
MJDP1 = MJDPGA + TOF1;

% N REV
N_rev = x(4);
N_rev2 = x(14);

% C3 launcher
v_inf_magn = x(6);
az = x(7);
el = x(8);

% GA stuff IN
v_inf_magn2 = x(9);
az2 = x(10);
el2 = x(11);
% GA stuff out
az3 = x(12);
el3 = x(13);

% chosing which asteroid to visit
IDP = x(5); %index of permutation, the column of the Permutation Matrix of the asteroids
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

% GA
MJDPGA_dim = MJDPGA*sim.TU/(3600*24);
[kep_GA,ksun] = uplanet(MJDPGA_dim, sim.ID_FLYBY);
[r_GA, v_GA] = sv_from_coe(kep_GA,ksun);
r_GA = r_GA/sim.DU;
v_GA = v_GA/sim.DU*sim.TU;

% passage at 1st ast
MJDP1_dim = MJDP1*sim.TU/(3600*24);
[kep_ast_1] = uNEO2(MJDP1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1, v1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
r1 = r1/sim.DU;
v1 = v1/sim.DU*sim.TU;

v_launcher = v_inf_magn*[cos(el)*cos(az); cos(el)*sin(az); sin(el)];
v_dep = v_EA + v_launcher;  %since parabolic escape (vinf = 0)

v_arr_GA = v_GA + v_inf_magn2*[cos(el2)*cos(az2); cos(el2)*sin(az2); sin(el2)];

[output_1] = NL_interpolator( r_EA , r_GA , v_dep , v_arr_GA , N_rev , TOFGA , sim.M ,sim.PS.Isp , sim );

v_dep_GA = v_GA + v_inf_magn2*[cos(el3)*cos(az3); cos(el3)*sin(az3); sin(el3)];

if sim.ID_FLYBY == 3
    RPlanet_flyby = astroConstants(23); % Radius_Earth, km
    muPlanet_flyby = astroConstants(13); % muEarth, km^3/s^2
    R_lim_from_planet = 500; % km, for earth is ok to avoid atmosphere
elseif sim.ID_FLYBY == 4
    RPlanet_flyby = astroConstants(24); % Radius_Earth, km
    muPlanet_flyby = astroConstants(14); % muEarth, km^3/s^2
    R_lim_from_planet = 200; % km, for mars is ok to avoid atmosphere
end
MJDPGA_dim = MJDPGA*sim.TU/(3600*24);
output.delta_v_p = flyby(RPlanet_flyby, muPlanet_flyby,R_lim_from_planet, ...
                         MJDPGA_dim, v_arr_GA, v_dep_GA, sim.ID_FLYBY);

if strcmp(string(output.delta_v_p), 'Not found')
    disp('Flyby  didn t converged')
end

M_after_GA = output_1.m(end);
[output_2] = NL_interpolator( r_GA , r1 , v_dep_GA , v1 , N_rev2 , TOF1 , M_after_GA ,sim.PS.Isp , sim );

output.mass_fract = (output_1.m(1) - output_2.m(end))/output_1.m(1);

output.dV_GA = sqrt((v_dep_GA(1) - v_arr_GA(1))^2+(v_dep_GA(2) - v_arr_GA(2))^2+(v_dep_GA(3) - v_arr_GA(3))^2);

%% Output encounter states
r_encounter.EA = r_EA;
r_encounter.GA = r_GA;
r_encounter.ast1 = r1;

v_encounter.EA = v_EA;
v_encounter.GA = v_GA;
v_encounter.ast1 = v1;

%% Output
output.t            = [output_1.t; output_1.t(end)+output_2.t];
output.m            = [output_1.m; output_2.m];
output.Thrust       = [output_1.Thrust; output_2.Thrust];
output.a            = [output_1.a; output_2.a];
output.r.leg1       = output_1.r; 
output.r.leg2       = output_2.r;
output.theta.leg1   = output_1.theta; 
output.theta.leg2   = output_2.theta;
output.z.leg1       = output_1.z;
output.z.leg2       = output_2.z;
output.Href.leg1    = output_1.Href;
output.Href.leg2    = output_2.Href;

end

