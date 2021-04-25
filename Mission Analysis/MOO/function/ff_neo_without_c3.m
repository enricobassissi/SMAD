function obj_fun = ff_neo_without_c3(x)
MJD01 = x(1);
MJDF1 = MJD01 + x(2);
buffer_time = x(3);
MJD02 = MJDF1 + buffer_time;
MJDF2 = MJD02 + x(4);

% Computing position and velocity of the planets in that days
% Departure from Earth
[kep_EA,ksun] = uplanet(MJD01, 3);
[r1, v1] = sv_from_coe(kep_EA,ksun);

% arrival at 1st ast
[r2,v2] = uNEO(MJDF1, '2020UE');

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
obj_fun(2) = TOF1 + TOF2;

end

