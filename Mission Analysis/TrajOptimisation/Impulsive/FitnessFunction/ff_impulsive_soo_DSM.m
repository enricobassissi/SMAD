function obj_fun = ff_impulsive_soo_DSM(x, data, sim)
% setting the input times
MJD01 = x(1);
TOF1a = x(10); % It's the time from Earth to DSM
TOF1b = x(2); % It's the time from DSM to asteroid 1
MJDF1 = MJD01 + TOF1a + TOF1b;
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

% DSM
rD_mag = x(11);
thetaD = x(12);

% Computing position and velocity of the planets in that days
% Departure from Earth
[kep_EA,ksun] = uplanet(MJD01, 3);
[r1, v1] = sv_from_coe(kep_EA,ksun); % km, km/s

% DSM
iEA = kep_EA(3);
rD = rD_mag*[cos(thetaD)*cos(iEA); sin(thetaD)*cos(iEA); sin(iEA)];

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
t1a_sec = MJD01*60*60*24;
t1b_sec = (MJD01 + TOF1a)*60*60*24;
t2_sec = MJDF1*60*60*24;
t3_sec = MJD02*60*60*24;
t4_sec = MJDF2*60*60*24;
t5_sec = MJD03*60*60*24;
t6_sec = MJDF3*60*60*24;
t7_sec = MJD04*60*60*24;
t8_sec = MJDF4*60*60*24;

% DV calculation with lambert

% Earth -> DSM
[~,~,~,~,VDd,VDa,~,~] = lambertMR(r1,rD,(t1b_sec - t1a_sec),ksun,0,0,0,0);

dv1 = sqrt((VDd(1)- v1(1))^2+(VDd(2)-v1(2))^2+(VDd(3)- v1(3))^2);

if dv1< sqrt(sim.C3_max) % vinf that the launcher can give max 
    dv_extra_launch = 0;
else
    c_launcher = 40; % penalty factor for dv_extra_launch
    dv_extra_launch = c_launcher*(dv1 - sqrt(sim.C3_max))^2; % penalty like, but not discard a priori
end

% DSM -> asteroid 1
[~,~,~,~,V1d,V1a,~,~] = lambertMR(rD,r2,(t2_sec - t1b_sec),ksun,0,0,0,0);
dv2 = sqrt((V1d(1)-VDa(1))^2+(V1d(2)-VDa(2))^2+(V1d(3)-VDa(3))^2);
dv3 = sqrt((V1a(1)- v2(1))^2+(V1a(2)-v2(2))^2+(V1a(3)-v2(3))^2);

if dv2 > 1
    CHECK_TERM_DSM =  50 + dv2^2;
else
    CHECK_TERM_DSM =  0;
end

[dv_tot45]=lambert_solver_rendezvous(r3,r4,v3,v4,t3_sec,t4_sec,ksun); 
[dv_tot67]=lambert_solver_rendezvous(r5,r6,v5,v6,t5_sec,t6_sec,ksun); 
[dv_tot89]=lambert_solver_rendezvous(r7,r8,v7,v8,t7_sec,t8_sec,ksun);



obj_fun = dv_extra_launch + dv2 + dv3 +dv_tot45 + dv_tot67 + dv_tot89 + CHECK_TERM_DSM;


end


