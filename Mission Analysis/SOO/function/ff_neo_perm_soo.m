function obj_fun = ff_neo_perm_soo(x, data, sim)
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
[kep_ast_1_1] = uNEO2(MJDF1,asteroid_1,data);
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

% % build of launcher help the injection, 3D;
% v_inf = v_inf_magn.*[cos(alpha)*cos(beta); 
%                     sin(alpha)*cos(beta); 
%                     sin(beta)];
% v_inj = v1+v_inf'; % v1 is row, v_inf would be column as built

% DV calculation with lambert
vlim = 100; %km/s
[~,dv1_12,dv2_12]=lambert_solver_rendezvous(r1,r2,v1,v2,t1_sec,t2_sec,ksun,vlim);
if dv1_12 < sqrt(sim.C3_max) % vinf that the launcher can give max
    [dv_tot34]=lambert_solver_rendezvous(r3,r4,v3,v4,t3_sec,t4_sec,ksun,vlim); 
    [dv_tot56]=lambert_solver_rendezvous(r5,r6,v5,v6,t5_sec,t6_sec,ksun,vlim); 
    [dv_tot78]=lambert_solver_rendezvous(r7,r8,v7,v8,t7_sec,t8_sec,ksun,vlim);

    obj_fun = dv2_12 + dv_tot34 + dv_tot56 + dv_tot78;
else
    obj_fun = NaN;
end

% TOF_tot = TOF1 + buffer_time1 + TOF2 + buffer_time2 + TOF3 + buffer_time3 + TOF4;
% 
% TOF_lim = 12*365; % 12 years
% if TOF_tot > TOF_lim
%     obj_fun = NaN;
% else

end

