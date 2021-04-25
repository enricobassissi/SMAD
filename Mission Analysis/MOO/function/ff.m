function obj_fun = ff(x)
MJD01 = x(1);
MJDF1 = MJD01 + x(2);
v_inf_magn = x(3);
alpha = x(4);
beta = x(5);
buffer_time = x(6);
MJD02 = MJDF1 + buffer_time;
MJDF2 = MJD02 + x(7);

% Computing position and velocity of the planets in that days
% Departure from Earth
[kep_EA,ksun] = uplanet(MJD01, 3);
[r1, v1] = sv_from_coe(kep_EA,ksun);

% arrival at mars
[kep_MA1,ksun] = uplanet(MJDF1, 4);
[r2, v2] = sv_from_coe(kep_MA1, ksun); 

% dep at mars after buffer
[kep_MA2,ksun] = uplanet(MJD02, 4);
[r3, v3] = sv_from_coe(kep_MA2, ksun); 

% arrival at jupiter
[kep_JU,ksun] = uplanet(MJDF2, 5);
[r4, v4] = sv_from_coe(kep_JU, ksun); 

% Converting mJ2000 in seconds
t1_sec = MJD01.*60.*60.*24;
t2_sec = MJDF1.*60.*60.*24;
t3_sec = MJD02.*60.*60.*24;
t4_sec = MJDF2.*60.*60.*24;

% build of launcher help the injection, 3D;
v_inf = v_inf_magn.*[cos(alpha)*cos(beta); 
                    sin(alpha)*cos(beta); 
                    sin(beta)];
v_inj = v1+v_inf'; % v1 is row, v_inf would be column as built

% DV calculation with lambert
[dv_tot1,~,~,TOF1]=lambert_solver_rendezvous( r1, r2, v_inj, v2, t1_sec, t2_sec, ksun, 50); % dvtot, tof

[dv_tot2,~,~,TOF2]=lambert_solver_rendezvous( r3, r4, v3, v4, t3_sec, t4_sec, ksun, 50); % dvtot, tof

obj_fun(1) = dv_tot1 + dv_tot2;
obj_fun(2) = TOF1 + TOF2;

end

