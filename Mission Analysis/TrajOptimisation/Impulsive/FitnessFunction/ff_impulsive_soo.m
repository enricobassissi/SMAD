function obj_fun = ff_impulsive_soo(x, data, sim)
% setting the input times
MJD01 = x(1);
TOF1 = x(2);
MJDF1 = MJD01 + TOF1;
buffer_time1 = x(3);
MJD02 = MJDF1 + buffer_time1;
TOF2 = x(4);
MJDF2 = MJD02 + TOF2;
buffer_time2 = x(6);
MJD03 = MJDF2 + buffer_time2;
TOF3 = x(7);
MJDF3 = MJD03 + TOF3; 
buffer_time3 = x(8);
MJD04 = MJDF3 + buffer_time3;
TOF4 =  x(9);
MJDF4 = MJD04 + TOF4; 

% chosing which asteroid to visit
IDP = round(x(5)); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_2 = data.PermutationMatrix(IDP,2);
asteroid_3 = data.PermutationMatrix(IDP,3);
asteroid_4 = data.PermutationMatrix(IDP,4);

% Computing position and velocity of the planets in that days
% Departure from Earth
[kep_EA,ksun] = uplanet(MJD01, 3);
[r1, v1] = sv_from_coe(kep_EA,ksun); % km, km/s
% arrival at 1st ast
[kep_ast_1_1] = uNEO2(MJDF1,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r2, v2] = sv_from_coe(kep_ast_1_1,ksun); % km, km/s

% dep at 1st after buffer
[kep_ast_1_2] = uNEO2(MJD02,asteroid_1,data);
[r3, v3] = sv_from_coe(kep_ast_1_2,ksun); % km, km/s
% arrival at 2nd asteroid 
[kep_ast_2_1] = uNEO2(MJDF2,asteroid_2,data);
[r4, v4] = sv_from_coe(kep_ast_2_1,ksun); % km, km/s

% dep at 2nd asteroid 
[kep_ast_2_2] = uNEO2(MJD03,asteroid_2,data);
[r5, v5] = sv_from_coe(kep_ast_2_2,ksun); % km, km/s
% arrival at 3rd asteroid 
[kep_ast_3_1] = uNEO2(MJDF3,asteroid_3,data);
[r6, v6] = sv_from_coe(kep_ast_3_1,ksun); % km, km/s

% dep at 3rd asteroid 
[kep_ast_3_2] = uNEO2(MJD04,asteroid_3,data);
[r7, v7] = sv_from_coe(kep_ast_3_2,ksun); % km, km/s
% arrival at 4th asteroid 
[kep_ast_4_1] = uNEO2(MJDF4,asteroid_4,data);
[r8, v8] = sv_from_coe(kep_ast_4_1,ksun); % km, km/s

% Converting mJ2000 in seconds
t1_sec = MJD01*60*60*24;
t2_sec = MJDF1*60*60*24;
t3_sec = MJD02*60*60*24;
t4_sec = MJDF2*60*60*24;
t5_sec = MJD03*60*60*24;
t6_sec = MJDF3*60*60*24;
t7_sec = MJD04*60*60*24;
t8_sec = MJDF4*60*60*24;

% DV calculation with lambert
[~,dv1_12,dv2_12]=lambert_solver_rendezvous(r1,r2,v1,v2,t1_sec,t2_sec,ksun);
if dv1_12 < sqrt(sim.C3_max) % vinf that the launcher can give max 
    dv_extra_launch = 0;
else
%     dv_extra_launch = dv1_12 - sqrt(sim.C3_max); %actually la differenza la paghi
    c_launcher = 30; % penalty factor for dv_extra_launch
    dv_extra_launch = c_launcher*(dv1_12 - sqrt(sim.C3_max))^2; % penalty like, but not discard a priori
end

[dv_tot34]=lambert_solver_rendezvous(r3,r4,v3,v4,t3_sec,t4_sec,ksun); 
[dv_tot56]=lambert_solver_rendezvous(r5,r6,v5,v6,t5_sec,t6_sec,ksun); 
[dv_tot78]=lambert_solver_rendezvous(r7,r8,v7,v8,t7_sec,t8_sec,ksun);

% Check of feasibility
CHECK_TERM_A = 0; CHECK_TERM_B = 0;
CHECK_TERM_C = 0; CHECK_TERM_D = 0;
tot_TOF = TOF1+TOF2+TOF3+TOF4;
if tot_TOF > 12*365
    CHECK_TERM_A = 100;
end
if dv_tot34 > 5
%     CHECK_TERM_B = 50;
    CHECK_TERM_B = (dv_tot34 - 6)^2;
end
if dv_tot56 > 5
%     CHECK_TERM_C = 50;
    CHECK_TERM_C = (dv_tot56 - 6)^2;
end
if dv_tot78 > 5
%     CHECK_TERM_D = 50;
    CHECK_TERM_D = (dv_tot78 - 6)^2;
end

CHECK_TERM = CHECK_TERM_A+CHECK_TERM_B+CHECK_TERM_C+CHECK_TERM_D;
% 
% pen_212 = 2;
% pen_34 = 1.5;
% pen_56 = 1;
% pen_78 = 0.8;
% obj_fun = dv_extra_launch + pen_212*dv2_12 + pen_34*dv_tot34 + ...
%     pen_56*dv_tot56 + pen_78*dv_tot78 + CHECK_TERM;

obj_fun = dv_extra_launch + dv2_12 + dv_tot34 + ...
    dv_tot56 + dv_tot78 + CHECK_TERM;

end

