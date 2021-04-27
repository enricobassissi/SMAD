function obj_fun = ff_neo_perm_adim(x, PermutationMatrix, sim)

% Measurament units definition
% sim.DU     Distance unit lenght[km]
% sim.TU     Time unit duration [s]

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
MJD01_dim = MJD01*sim.TU/(3600*24); % days
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r1, v1] = sv_from_coe(kep_EA,ksun);
r1 = r1/sim.DU; v1 = v1/sim.DU*sim.TU;
% arrival at 1st ast
MJDF1_dim = MJDF1*sim.TU/(3600*24);
[r2,v2] = uNEO(MJDF1_dim, asteroid_1);
r2 = r2/sim.DU; v2 = v2/sim.DU*sim.TU;

% dep at 1st after buffer
MJD02_dim = MJD02*sim.TU/(3600*24);
[r3, v3] = uNEO(MJD02_dim, asteroid_1);
r3 = r3/sim.DU; v3 = v3/sim.DU*sim.TU;
% arrival at 2nd asteroid 
MJDF2_dim = MJDF2*sim.TU/(3600*24);
[r4, v4] = uNEO(MJDF2_dim, asteroid_2);
r4 = r4/sim.DU; v4 = v4/sim.DU*sim.TU;

% dep at 2nd asteroid 
MJD03_dim = MJD03*sim.TU/(3600*24);
[r5, v5] = uNEO(MJD03_dim, asteroid_2);
r5 = r5/sim.DU; v5 = v5/sim.DU*sim.TU;
% arrival at 3rd asteroid 
MJDF3_dim = MJDF3*sim.TU/(3600*24);
[r6, v6] = uNEO(MJDF3_dim, asteroid_3);
r6 = r6/sim.DU; v6 = v6/sim.DU*sim.TU;

% dep at 3rd asteroid 
MJD04_dim = MJD04*sim.TU/(3600*24);
[r7, v7] = uNEO(MJD04_dim, asteroid_3);
r7 = r7/sim.DU; v7 = v7/sim.DU*sim.TU;
% arrival at 4th asteroid 
MJDF4_dim = MJDF4*sim.TU/(3600*24);
[r8, v8] = uNEO(MJDF4_dim, asteroid_4);
r8 = r8/sim.DU; v8 = v8/sim.DU*sim.TU;

% build of launcher help the injection, 3D;
v_inf = v_inf_magn.*[cos(alpha)*cos(beta); 
                    sin(alpha)*cos(beta); 
                    sin(beta)];
v_inj = v1+v_inf'; % v1 is row, v_inf would be column as built

% DV calculation with lambert
vlim = 100;
[dv_tot1]=lambert_solver_rendezvous(r1,r2,v_inj,v2,MJD01,MJDF1,sim.mu,vlim);
[dv_tot2]=lambert_solver_rendezvous(r3,r4,v3,v4,MJD02,MJDF2,sim.mu,vlim); 
[dv_tot3]=lambert_solver_rendezvous(r5,r6,v5,v6,MJD03,MJDF3,sim.mu,vlim); 
[dv_tot4]=lambert_solver_rendezvous(r7,r8,v7,v8,MJD04,MJDF4,sim.mu,vlim);

obj_fun(1) = (dv_tot1 + dv_tot2 + dv_tot3 + dv_tot4)*sim.DU/sim.TU; % km/s
   
if isnan(obj_fun(1))
    obj_fun(2) = NaN;
else
    obj_fun(2) = (TOF1 + buffer_time1 + TOF2 + buffer_time2 + TOF3 + buffer_time3 + TOF4)*sim.TU/(3600*24); %days
end

end

