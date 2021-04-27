function obj_fun = ff_neo_perm_soo(x,PermutationMatrix)
% setting the input times
MJD01 = x(1);
TOF1 = x(2);
MJDF1 = MJD01 + TOF1;
buffer_time1 = x(6);
MJD02 = MJDF1 + buffer_time1;
TOF2 = x(7);
MJDF2 = MJD02 + TOF2;
buffer_time2 = x(9);
MJD03 = MJDF2 + buffer_time2;
TOF3 = x(10);
MJDF3 = MJD03 + TOF3; 
buffer_time3 = x(11);
MJD04 = MJDF3 + buffer_time3;
TOF4 =  x(12);
MJDF4 = MJD04 + TOF4; 

% launcher departure vector, that set the required c3
v_inf_magn = x(3);
alpha = x(4);
beta = x(5);

% chosing which asteroid to visit
IDP = round(x(8)); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = PermutationMatrix(IDP,1);
asteroid_2 = PermutationMatrix(IDP,2);
asteroid_3 = PermutationMatrix(IDP,3);
asteroid_4 = PermutationMatrix(IDP,4);

% Computing position and velocity of the planets in that days
% Departure from Earth
[kep_EA,ksun] = uplanet(MJD01, 3);
[r1, v1] = sv_from_coe(kep_EA,ksun);
% arrival at 1st ast
[r2,v2] = uNEO(MJDF1, asteroid_1);

% dep at 1st after buffer
[r3, v3] = uNEO(MJD02, asteroid_1);
% arrival at 2nd asteroid 
[r4, v4] = uNEO(MJDF2, asteroid_2);

% dep at 2nd asteroid 
[r5, v5] = uNEO(MJD03, asteroid_2);
% arrival at 3rd asteroid 
[r6, v6] = uNEO(MJDF3, asteroid_3);

% dep at 3rd asteroid 
[r7, v7] = uNEO(MJD04, asteroid_3);
% arrival at 4th asteroid 
[r8, v8] = uNEO(MJDF4, asteroid_4);

% Converting mJ2000 in seconds
t1_sec = MJD01*60*60*24;
t2_sec = MJDF1*60*60*24;
t3_sec = MJD02*60*60*24;
t4_sec = MJDF2*60*60*24;
t5_sec = MJD03*60*60*24;
t6_sec = MJDF3*60*60*24;
t7_sec = MJD04*60*60*24;
t8_sec = MJDF4*60*60*24;

% build of launcher help the injection, 3D;
v_inf = v_inf_magn.*[cos(alpha)*cos(beta); 
                    sin(alpha)*cos(beta); 
                    sin(beta)];
v_inj = v1+v_inf'; % v1 is row, v_inf would be column as built

% DV calculation with lambert
vlim = 100;
[dv_tot1]=lambert_solver_rendezvous(r1,r2,v_inj,v2,t1_sec,t2_sec,ksun,vlim);
[dv_tot2]=lambert_solver_rendezvous(r3,r4,v3,v4,t3_sec,t4_sec,ksun,vlim); 
[dv_tot3]=lambert_solver_rendezvous(r5,r6,v5,v6,t5_sec,t6_sec,ksun,vlim); 
[dv_tot4]=lambert_solver_rendezvous(r7,r8,v7,v8,t7_sec,t8_sec,ksun,vlim);

obj_fun(1) = dv_tot1 + dv_tot2 + dv_tot3 + dv_tot4;

TOF_tot = TOF1 + buffer_time1 + TOF2 + buffer_time2 + TOF3 + buffer_time3 + TOF4;

TOF_lim = 12*365; % 12 years
if TOF_tot > TOF_lim
    obj_fun = NaN;
else

end

