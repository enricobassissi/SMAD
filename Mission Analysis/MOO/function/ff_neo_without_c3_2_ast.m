function obj_fun = ff_neo_without_c3_2_ast(x)
MJD01 = x(1);
MJDF1 = MJD01 + x(2);

% Computing position and velocity of the planets in that days
% Departure from Earth
[~,ksun] = uplanet(MJD01, 3);
% [r1, v1] = sv_from_coe(kep_EA,ksun); 
[r1, v1] = uNEO(MJD01, '2006HX57');

% arrival at 1st ast
[r2,v2] = uNEO(MJDF1, '2020UE');

% Converting mJ2000 in seconds
t1_sec = MJD01.*60.*60.*24;
t2_sec = MJDF1.*60.*60.*24;

% DV calculation with lambert
vlim=1e2;
[dv_tot1,~,~,TOF1]=lambert_solver_rendezvous( r1, r2, v1, v2, t1_sec, t2_sec, ksun, vlim); % dvtot, tof

obj_fun(1) = dv_tot1 ;
obj_fun(2) = TOF1 ;

end

