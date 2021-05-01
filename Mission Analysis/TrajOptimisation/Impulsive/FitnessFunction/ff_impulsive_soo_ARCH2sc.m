function obj_fun = ff_impulsive_soo_ARCH2sc(x, data, sim)

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

% 1ST SPACECRAFT ASTEROID OBJECTIVES
% choosing which asteroid to visit
IDP1 = round(x(6)); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP1,1);
asteroid_2 = data.PermutationMatrix(IDP1,2);

% 2ND SPACECRAFT ASTEROID OBJECTIVES
IDP2 = round(x(7)); % index for 2nd permutation matrix to be built inside depending on the first 2 selected asteroids
asteroid_sequence = [asteroid_1,asteroid_2];
TF = contains(data.asteroid_names,asteroid_sequence);
not_asteroid_sequence = data.asteroid_names(~TF);
% HowMany_for2ndSC = factorial(length(not_asteroid_sequence)) / factorial(length(not_asteroid_sequence) - 2);
[PermutationMatrix_SC2, ~] = permnUnique(not_asteroid_sequence, 2);
asteroid_a = PermutationMatrix_SC2(IDP2,1);
asteroid_b = PermutationMatrix_SC2(IDP2,2);

% Computing position and velocity of the planets in that days
% Departure from Earth
[kep_EA,ksun] = uplanet(MJD01, 3);
[rEA, vEA] = sv_from_coe(kep_EA,ksun); % km, km/s
% passage of 1st sc at 1st ast
[kep_ast_1] = uNEO2(MJDP1,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1, v1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
% passage of 1st sc at 2nd asteroid 
[kep_ast_2] = uNEO2(MJDP2,asteroid_2,data);
[r2, v2] = sv_from_coe(kep_ast_2,ksun); % km, km/s

% passage of 2nd sc at 1st ast
[kep_ast_a] = uNEO2(MJDPa,asteroid_a,data);
[ra, va] = sv_from_coe(kep_ast_a,ksun); % km, km/s
% passage of 2nd sc at 2nd ast
[kep_ast_b] = uNEO2(MJDPb,asteroid_b,data);
[rb, vb] = sv_from_coe(kep_ast_b,ksun); % km, km/s

% Converting TOFs in seconds for Lambert
ToF_EAast1_sec = TOF1*60*60*24;
ToF_ast12_sec = TOF2*60*60*24;

ToF_EAasta_sec = TOFa*60*60*24;
ToF_astab_sec = TOFb*60*60*24;

% DV calculation with lambert
% SPACECRAFT 1
% Earth -> asteroid 1
[~,~,~,~,VI_EAast1,VF_EAast1,~,~] = lambertMR(rEA,r1,ToF_EAast1_sec,ksun,0,0,0,0);
dv1_EAast1 = sqrt((VI_EAast1(1)-vEA(1))^2+(VI_EAast1(2)-vEA(2))^2+(VI_EAast1(3)-vEA(3))^2);
if dv1_EAast1 < sqrt(sim.C3_max) % vinf that the launcher can give max 
    dv_extra_launch1st = 0;
else
    % actually you would pay the difference, so we put a very big number to
    % let the optimizer to not consider this solution because too expansive
    % dv1_EAast1 - sqrt(sim.C3_max); 
    dv_extra_launch1st = 30; % very high number, arbitrary
end
dv2_EAast1 = sqrt((VF_EAast1(1)-v1(1))^2+(VF_EAast1(2)-v1(2))^2+(VF_EAast1(3)-v1(3))^2);

% asteroid 1 -> 2
[~,~,~,~,VI_ast12,VF_ast12,~,~] = lambertMR(r1,r2,ToF_ast12_sec,ksun,0,0,0,0);
dv1_ast12 = sqrt((VI_ast12(1)-v1(1))^2+(VI_ast12(2)-v1(2))^2+(VI_ast12(3)-v1(3))^2);
dv2_ast12 = sqrt((VF_ast12(1)-v2(1))^2+(VF_ast12(2)-v2(2))^2+(VF_ast12(3)-v2(3))^2);

% dV of flyby passage on asteroid 1
dv_passage_ast1 = sqrt((VI_ast12(1)-VF_EAast1(1))^2+(VI_ast12(2)-VF_EAast1(2))^2+(VI_ast12(3)-VF_EAast1(3))^2);

% SPACECRAFT 2
% Earth -> asteroid a
[~,~,~,~,VI_EAasta,VF_EAasta,~,~] = lambertMR(rEA,ra,ToF_EAasta_sec,ksun,0,0,0,0);
dv1_EAasta = sqrt((VI_EAasta(1)-vEA(1))^2+(VI_EAasta(2)-vEA(2))^2+(VI_EAasta(3)-vEA(3))^2);
if dv1_EAasta < sqrt(sim.C3_max) % vinf that the launcher can give max 
    dv_extra_launch2nd = 0;
else
    % actually you would pay the difference, so we put a very big number to
    % let the optimizer to not consider this solution because too expansive
    % dv1_EAast1 - sqrt(sim.C3_max); 
    dv_extra_launch2nd = 30; % very high number, arbitrary
end
dv2_EAasta = sqrt((VF_EAasta(1)-va(1))^2+(VF_EAasta(2)-va(2))^2+(VF_EAasta(3)-va(3))^2);

% asteroid a -> b
[~,~,~,~,VI_astab,VF_astab,~,~] = lambertMR(ra,rb,ToF_astab_sec,ksun,0,0,0,0);
dv1_astab = sqrt((VI_astab(1)-va(1))^2+(VI_astab(2)-va(2))^2+(VI_astab(3)-va(3))^2);
dv2_astab = sqrt((VF_astab(1)-vb(1))^2+(VF_astab(2)-vb(2))^2+(VF_astab(3)-vb(3))^2);

% dV of flyby passage on asteroid a
dv_passage_asta = sqrt((VI_astab(1)-VF_EAasta(1))^2+(VI_astab(2)-VF_EAasta(2))^2+(VI_astab(3)-VF_EAasta(3))^2);

% if the last dv is a flyby it can go wherever it wants after the last encounter, both the sc
obj_fun = dv_extra_launch1st + dv_extra_launch2nd + dv_passage_ast1 + dv_passage_asta; 

% else if the last dv is a rendezvous with the last object
% obj_fun = all previous + 2 dv of rendezvous of the 2 sc

% obj_fun = TOF1+TOF2+TOF3+TOF4;

end

