function obj_fun = ff_LT_1asteroid(x,sim,data)
% setting the input times
MJD01 = x(1);
TOF1 = x(2);
MJDF1 = MJD01 + TOF1;

% Computing position and velocity of the planets in that days
% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r1, v1] = sv_from_coe(kep_EA,ksun);

r1 = r1/sim.DU;
v1 = v1/sim.DU*sim.TU;

% passage at 1st ast
MJDF1_dim = MJDF1*sim.TU/(3600*24);
[kep_ast_1] = uNEO2(MJDF1_dim,sim.ast1,data); % [km,-,rad,rad,rad,wrapped rad]
[r2, v2] = sv_from_coe(kep_ast_1,ksun); % km, km/s
% adimensionalise
r2 = r2/sim.DU;
v2 = v2/sim.DU*sim.TU;

% N REV
N_rev= round(x(3));


% vinf
vinf_mag = x(4);
alpha = x(5);
beta = x(6);

vdep = v1 + vinf_mag*[cos(alpha)*cos(beta); sin(alpha)*cos(beta); sin(beta)] ;  % ok?
%vdep = v1;
varr = v2; %since we want to rendez-vous

[output] = NL_interpolator( r1 , r2 , vdep , varr , N_rev , TOF1 ,sim.M ,sim.PS.Is ,sim );

%obj_fun(1) = TOF1;
   
obj_fun = (output.m(1) - output.m(end))/output.m(1); % la massa è dimenionale

end