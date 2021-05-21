function obj_fun = ff_ea_1ast_LT_RV_GA_soo_NLI(x,sim,data)
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

% initialise and if everything fine it rewrite it well
obj_fun = 10;

if ~isnan(output_1.Thrust) % if is not nan
    
    v_dep_GA = v_GA + v_inf_magn2*[cos(el3)*cos(az3); cos(el3)*sin(az3); sin(el3)];
    
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
    MJDPGA_dim = MJDPGA*sim.TU/(3600*24);
    v_arr_GA_dim = v_arr_GA.*sim.DU./sim.TU;
    v_dep_GA_dim = v_dep_GA.*sim.DU./sim.TU;
    [delta_v_p,rp] = flyby(RPlanet_flyby, muPlanet_flyby,R_lim_from_planet, ...
                      MJDPGA_dim, v_arr_GA_dim, v_dep_GA_dim, sim.ID_FLYBY, R_SOI_PL);
    
    if ~strcmp(string(delta_v_p), 'Not found')
        
        M_after_GA = output_1.m(end);
        [output_2] = NL_interpolator( r_GA , r1 , v_dep_GA , v1 , N_rev2 , TOF1 , M_after_GA ,sim.PS.Isp , sim );

        if ~isnan(output_2.Thrust) % if is not nan

            mass_fract = (output_1.m(1) - output_2.m(end))/output_1.m(1);

            % put one after the other, all the thrust profiles
            T_append = [output_1.Thrust(:,1),output_1.Thrust(:,2),output_1.Thrust(:,3);
                        output_2.Thrust(:,1),output_2.Thrust(:,2),output_2.Thrust(:,3)]; 
            T = sqrt(T_append(:,1).^2 + T_append(:,3).^2);

            if abs(max(T)) <= 0.025 % bepi colombo is 250 mN
                if mass_fract > 0 && mass_fract < 1 
                    obj_fun = mass_fract;
                    disp('success')
                else
                    obj_fun = obj_fun + 10;
                end
            else
                obj_fun = obj_fun + 21;
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

end

