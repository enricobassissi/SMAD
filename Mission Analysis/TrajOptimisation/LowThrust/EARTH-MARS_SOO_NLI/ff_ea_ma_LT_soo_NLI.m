function obj_fun = ff_ea_ma_LT_soo_NLI(x,sim)
% setting the input times
MJD01 = x(1);
TOF1 = x(2);
%TOF1 = 1000;
MJDF1 = MJD01 + TOF1;

% Computing position and velocity of the planets in that days
% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r1, v1] = sv_from_coe(kep_EA,ksun);

r1 = r1/sim.DU;
v1 = v1/sim.DU*sim.TU;

% arrival on mars
MJDF1_dim = MJDF1*sim.TU/(3600*24);
[kep_MA,ksun] = uplanet(MJDF1_dim, 4);
[r2, v2] = sv_from_coe(kep_MA,ksun);
% adimensionalise
r2 = r2/sim.DU;
v2 = v2/sim.DU*sim.TU;

% N REV
N_rev= round(x(3));

vdep = v1;  %since parabolic escape (vinf = 0)
varr = v2; %since we want to rendez-vous

Isp = sim.PS.Isp;

[output] = NL_interpolator( r1 , r2 , vdep , varr , N_rev , TOF1 , sim.M ,Isp ,sim );

obj_fun = (output.m(1) - output.m(end))/output.m(1); % la massa è dimenionale

end

