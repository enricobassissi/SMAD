MJD01=9000;
MJDF1=9200;
MJD02=9350;
MJDF2=9750;

[kep_EA,ksun] = uplanet(MJD01, 3);
[r1, v1] = sv_from_coe(kep_EA,ksun); %km, km/s

% arrival at 1st ast
[r2,v2] = uNEO(MJDF1, '2020UE'); %km, km/s

% dep at 1st after buffer
[r3, v3] = uNEO(MJD02, '2020UE');

% arrival at 2nd asteroid 
[r4, v4] = uNEO(MJDF2, '2006HX57');

% Converting mJ2000 in seconds
t1_sec = MJD01.*60.*60.*24;
t2_sec = MJDF1.*60.*60.*24;
t3_sec = MJD02.*60.*60.*24;
t4_sec = MJDF2.*60.*60.*24;

% DV calculation with lambert
vlim=1e2;
[dv_tot1,~,~,TOF1]=lambert_solver_rendezvous( r1, r2, v1, v2, t1_sec, t2_sec, ksun, vlim); % dvtot, tof

[dv_tot2,~,~,TOF2]=lambert_solver_rendezvous( r3, r4, v3, v4, t3_sec, t4_sec, ksun, vlim); % dvtot, tof

obj_fun(1) = dv_tot1 + dv_tot2;
obj_fun(2) = TOF1 + TOF2