function obj_fun = ff_impulsive_soo_ARCH1plus4(x, data, sim)
% setting the input times
MJD01 = x(1);
TOF1 = x(2);
MJDP1 = MJD01 + TOF1; % mjd2000 passage on ast 1
TOF2 = x(3);
MJDP2 = MJDP1 + TOF2;
TOF3 = x(4);
MJDP3 = MJDP2 + TOF3; 
TOF4 =  x(5);
MJDP4 = MJDP3 + TOF4; 

% chosing which asteroid to visit
IDP = round(x(6)); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_2 = data.PermutationMatrix(IDP,2);
asteroid_3 = data.PermutationMatrix(IDP,3);
asteroid_4 = data.PermutationMatrix(IDP,4);

% Computing position and velocity of the planets in that days
% Departure from Earth
[kep_EA,ksun] = uplanet(MJD01, 3);
[rEA, vEA] = sv_from_coe(kep_EA,ksun); % km, km/s
% passage at 1st ast
[kep_ast_1] = uNEO2(MJDP1,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1, v1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
% passage at 2nd asteroid 
[kep_ast_2] = uNEO2(MJDP2,asteroid_2,data);
[r2, v2] = sv_from_coe(kep_ast_2,ksun); % km, km/s
% passage at 3rd asteroid 
[kep_ast_3] = uNEO2(MJDP3,asteroid_3,data);
[r3, v3] = sv_from_coe(kep_ast_3,ksun); % km, km/s
% arrival at 4th asteroid 
[kep_ast_4] = uNEO2(MJDP4,asteroid_4,data);
[r4, v4] = sv_from_coe(kep_ast_4,ksun); % km, km/s

% Converting mJ2000 in seconds
t1_sec = MJD01*60*60*24;
t2_sec = MJDP1*60*60*24;
t3_sec = MJDP2*60*60*24;
t4_sec = MJDP3*60*60*24;
t5_sec = MJDP4*60*60*24;

% DV calculation with lambert
[~,dv1_EA_1,dv2_EA_1]=lambert_solver_flyby(rEA,r1,vEA,v1,t1_sec,t2_sec,ksun);
if dv1_EA_1 < sqrt(sim.C3_max) % vinf that the launcher can give max 
    dv_extra_launch = 0;
else
    dv_extra_launch = dv1_EA_1 - sqrt(sim.C3_max); %actually la differenza la paghi
end

% asteroid 1 -> 2
[~,dv1_1_2,dv2_1_2]=lambert_solver_flyby(r1,r2,v1,v2,t2_sec,t3_sec,ksun); 
% asteroid 2 -> 3
[~,dv1_2_3,dv2_2_3]=lambert_solver_flyby(r2,r3,v2,v3,t3_sec,t4_sec,ksun); 
% asteroid 3 -> 4
[~,dv1_3_4,~]=lambert_solver_flyby(r3,r4,v3,v4,t4_sec,t5_sec,ksun);

obj_fun = dv_extra_launch + abs(dv2_EA_1-dv1_1_2) + abs(dv2_1_2-dv1_2_3) + abs(dv2_2_3-dv1_3_4); % because the last dv is a flyby it can go wherever it wants
% obj_fun = TOF1+TOF2+TOF3+TOF4;

end

