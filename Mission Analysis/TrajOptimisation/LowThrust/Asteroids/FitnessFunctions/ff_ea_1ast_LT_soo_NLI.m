function obj_fun = ff_ea_1ast_LT_soo_NLI(x,sim,data)
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
% varr = v2; %since we want to rendezvous

[output] = NL_interpolator( r_EA , r1 , v_dep , v1 , N_rev , TOF1 , sim.M ,sim.PS.Isp , sim );

obj_fun = 10;
tol_TOF = 1;
if ~isnan(output.Thrust(1,1)) && abs(output.t(end) - TOF1) < tol_TOF % if is not nan
    
    mass_fract = (output.m(1) - output.m(end))/output.m(1);
    
    T = sqrt(output.Thrust(:,1).^2 + output.Thrust(:,3).^2);
    
    if abs(max(T)) < sim.max_Available_Thrust % bepi colombo is 250 mN
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

end

