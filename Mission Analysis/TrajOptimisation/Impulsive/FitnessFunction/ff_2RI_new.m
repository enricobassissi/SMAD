function obj_fun = ff_2RI_new(x, data, sim)

% Nomenclature
% quantities for SpaceCraft1 will have numbers (1,2,...)
% quantities for SC 2 will have letters (a,b,...)

% setting the input times
MJD01 = x(1); % departure time for both the sc

% 1st spacecraft characteristic times
TOF1 = x(2); % tof sc1 to 1st asteroid
MJDF1 = MJD01 + TOF1; % mjd2000 arrival of 1st sc on ast 1
buffer_time1 = x(3); % time to spend on ast 1 by sc 1
MJD02 = MJDF1 + buffer_time1; % departure of sc 1 from ast 1
TOF2 = x(4); % tof sc1 to 2nd asteroid
MJDF2 = MJD02 + TOF2; % mjd2000 passage of 1st sc on ast 2

% 2nd spacecraft characteristic times
TOFa = x(5); % tof sc2 to 1st asteroid
MJDFa = MJD01 + TOFa; % arrival time of sc 2 on ast a
buffer_time2 = x(6); % time to spend on ast a by sc 2
MJD0a = MJDFa + buffer_time2; % departure of sc 2 from ast a
TOFb =  x(7); % tof sc2 to 2nd asteroid
MJDFb = MJD0a + TOFb; 

% 1ST SPACECRAFT ASTEROID OBJECTIVES
% choosing which asteroid to visit
IDP1 = round(x(8)); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP1,1);
asteroid_2 = data.PermutationMatrix(IDP1,2);

% 2ND SPACECRAFT ASTEROID OBJECTIVES
% IDP2 = round(x(9)); % index for 2nd permutation matrix to be built inside depending on the first 2 selected asteroids
IDP_temp_2 = round(x(9)); % index for 2nd permutation matrix to be built inside depending on the first 2 selected asteroids
asteroid_sequence = [asteroid_1,asteroid_2];
TF = contains(data.asteroid_names,asteroid_sequence);
data_elements_matrix_2SC = data.data_elements_matrix(~TF,:);
[~, PermutationMatrix_2SC, HowMany_2SC] = ...
            sequences_local_pruning2(data_elements_matrix_2SC, data.p_number);
IDP2 = ceil(IDP_temp_2*HowMany_2SC/100);
asteroid_a = PermutationMatrix_2SC(IDP2,1);
asteroid_b = PermutationMatrix_2SC(IDP2,2);

% Computing position and velocity of the planets in that days
% Departure from Earth
[kep_EA,ksun] = uplanet(MJD01, 3);
[rEA, vEA] = sv_from_coe(kep_EA,ksun); % km, km/s
% arrival of 1st sc at 1st ast
[kep_ast_11] = uNEO2(MJDF1,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1_arr, v1_arr] = sv_from_coe(kep_ast_11,ksun); % km, km/s
% departure of 1st sc at 1st ast
[kep_ast_12] = uNEO2(MJD02,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1_dep, v1_dep] = sv_from_coe(kep_ast_12,ksun); % km, km/s
% arrival of 1st sc at 2nd asteroid 
[kep_ast_2] = uNEO2(MJDF2,asteroid_2,data);
[r2, v2] = sv_from_coe(kep_ast_2,ksun); % km, km/s

% arrival of 2nd sc at 1st ast
[kep_ast_a1] = uNEO2(MJDFa,asteroid_a,data);
[ra_arr, va_arr] = sv_from_coe(kep_ast_a1,ksun); % km, km/s
% departure of 2nd sc at 1st ast
[kep_ast_a2] = uNEO2(MJD0a,asteroid_a,data);
[ra_dep, va_dep] = sv_from_coe(kep_ast_a2,ksun); % km, km/s
% arrival of 2nd sc at 2nd ast
[kep_ast_b] = uNEO2(MJDFb,asteroid_b,data);
[rb, vb] = sv_from_coe(kep_ast_b,ksun); % km, km/s

% Converting TOFs in seconds for Lambert
ToF_EAast1_sec = TOF1*60*60*24;
buffer_time1_sec = buffer_time1*60*60*24;
ToF_ast12_sec = TOF2*60*60*24;

ToF_EAasta_sec = TOFa*60*60*24;
buffer_time2_sec = buffer_time2*60*60*24;
ToF_astab_sec = TOFb*60*60*24;

% DV calculation with lambert
% SPACECRAFT 1
% Earth -> asteroid 1
[~,~,~,~,VI_EAast1,VF_EAast1,~,~] = lambertMR(rEA,r1_arr,ToF_EAast1_sec,ksun,0,0,0,0);
dv1_EAast1 = sqrt((VI_EAast1(1)-vEA(1))^2+(VI_EAast1(2)-vEA(2))^2+(VI_EAast1(3)-vEA(3))^2);
if dv1_EAast1 < sqrt(sim.C3_max) % vinf that the launcher can give max 
    dv_extra_launch1st = 0;
else
    % actually you would pay the difference, so we put a very big number to
    % let the optimizer to not consider this solution because too expansive
    % dv1_EAast1 - sqrt(sim.C3_max); 
%     dv_extra_launch1st = 30; % very high number, arbitrary
    dv_extra_launch1st = 10*(dv1_EAast1 - sqrt(sim.C3_max))^2;
end
dv2_EAast1 = sqrt((VF_EAast1(1)-v1_arr(1))^2+(VF_EAast1(2)-v1_arr(2))^2+(VF_EAast1(3)-v1_arr(3))^2);

% asteroid 1 -> 2
[~,~,~,~,VI_ast12,VF_ast12,~,~] = lambertMR(r1_dep,r2,ToF_ast12_sec,ksun,0,0,0,0);
dv1_ast12 = sqrt((VI_ast12(1)-v1_dep(1))^2+(VI_ast12(2)-v1_dep(2))^2+(VI_ast12(3)-v1_dep(3))^2);
dv2_ast12 = sqrt((VF_ast12(1)-v2(1))^2+(VF_ast12(2)-v2(2))^2+(VF_ast12(3)-v2(3))^2);

% dV total of rendezvous of sc 1: entering ast 1 orbit, exiting 1, entering 2
dv_rv_sc1 = dv2_EAast1 + dv1_ast12 + dv2_ast12;

% SPACECRAFT 2
% Earth -> asteroid a
[~,~,~,~,VI_EAasta,VF_EAasta,~,~] = lambertMR(rEA,ra_arr,ToF_EAasta_sec,ksun,0,0,0,0);
dv1_EAasta = sqrt((VI_EAasta(1)-vEA(1))^2+(VI_EAasta(2)-vEA(2))^2+(VI_EAasta(3)-vEA(3))^2);
if dv1_EAasta < sqrt(sim.C3_max) % vinf that the launcher can give max 
    dv_extra_launch2nd = 0;
else
    % actually you would pay the difference, so we put a very big number to
    % let the optimizer to not consider this solution because too expansive
    % dv1_EAast1 - sqrt(sim.C3_max); 
%     dv_extra_launch2nd = 30; % very high number, arbitrary
    dv_extra_launch2nd = 10*(dv1_EAasta - sqrt(sim.C3_max))^2;
end
dv2_EAasta = sqrt((VF_EAasta(1)-va_arr(1))^2+(VF_EAasta(2)-va_arr(2))^2+(VF_EAasta(3)-va_arr(3))^2);

% asteroid a -> b
[~,~,~,~,VI_astab,VF_astab,~,~] = lambertMR(ra_dep,rb,ToF_astab_sec,ksun,0,0,0,0);
dv1_astab = sqrt((VI_astab(1)-va_dep(1))^2+(VI_astab(2)-va_dep(2))^2+(VI_astab(3)-va_dep(3))^2);
dv2_astab = sqrt((VF_astab(1)-vb(1))^2+(VF_astab(2)-vb(2))^2+(VF_astab(3)-vb(3))^2);

% dV total of rendezvous of sc 2: entering ast a orbit, exiting a, entering b
dv_rv_sc2 = dv2_EAasta + dv1_astab + dv2_astab;

% if the last dv is a flyby it can go wherever it wants after the last encounter, both the sc
% obj_fun = dv_extra_launch1st + dv_extra_launch2nd + dv_passage_ast1 + dv_passage_asta; 

% else if the last dv is a rendezvous with the last object
% obj_fun = all previous + 2 dv of rendezvous of the 2 sc

% Check of feasibility
CHECK_TERM = 0; CHECK_TERM_TOF = 0;
CHECK_TERM_A = 0; CHECK_TERM_B = 0;
CHECK_TERM_C = 0; CHECK_TERM_D = 0;
tot_TOF = TOF1+TOF2+TOFa+TOFb;
if tot_TOF > 12*365
    CHECK_TERM_TOF = 100;
end
if dv2_EAast1 > 7
    CHECK_TERM_A = (dv2_EAast1)^2;
end
if dv2_ast12 > 7
    CHECK_TERM_B = (dv2_ast12)^2;
end
if dv2_EAasta > 7
    CHECK_TERM_C = (dv2_EAasta)^2;
end
if dv2_astab > 7
    CHECK_TERM_D = (dv2_astab)^2;
end
CHECK_TERM = max([CHECK_TERM_A,CHECK_TERM_B,CHECK_TERM_C,CHECK_TERM_D]);

% take into account also the rel velocity at the asteroids
% c = 0.01; % penalty function
% weight less the vrel
% avg_dVrel_passage = mean([dv2_EAast1+dv2_ast12+dv2_EAasta+dv2_astab]);
% obj_fun = dv_extra_launch1st + dv_extra_launch2nd + dv_passage_ast1 + dv_passage_asta +...
%     c*((dv2_EAast1-avg_dVrel_passage)^2 + (dv2_ast12-avg_dVrel_passage)^2 +...
%     (dv2_EAasta-avg_dVrel_passage)^2 + (dv2_astab-avg_dVrel_passage)^2) + CHECK_TERM;

% penalty related to different mass of two spacecraft
c1 = 0.05; % penalty function
penDELTAdv = c1*(dv_rv_sc1- dv_rv_sc2)^2;
obj_fun = dv_extra_launch1st + dv_extra_launch2nd + dv_rv_sc1 + dv_rv_sc2 +...
    penDELTAdv + CHECK_TERM_TOF + 0.5*CHECK_TERM; 

end

