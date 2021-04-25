function obj_fun = ff_EA2MA_2(x)
MJD0 = x(1);
MJDF = x(1)+x(2);
v_inf_magn = x(3);
alpha = x(4);
beta = x(5);

% Computing position and velocity of the planets in that days
[kep_EA,ksun] = uplanet(MJD0, 3);
[r1, v1] = sv_from_coe(kep_EA,ksun);

[kep_MA,ksun] = uplanet(MJDF, 4);
[r2, v2] = sv_from_coe(kep_MA, ksun); 


% Converting mJ2000 in seconds
t1_sec=MJD0.*60.*60.*24;
t2_sec=MJDF.*60.*60.*24;

% build of launcher help the injection, 3D;
v_inf = v_inf_magn.*[cos(alpha)*cos(beta); 
                    sin(alpha)*cos(beta); 
                    sin(beta)];
v_inj = v1+v_inf'; % v1 is row, v_inf would be column as built

% DV calculation with Lambert, and tof
[obj_fun(1),~,~,obj_fun(2)]=lambert_solver_rendezvous( r1, r2, v_inj, v2, t1_sec, t2_sec, ksun, 50); % dvtot, tof

end

