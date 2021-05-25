function obj_fun = ff_ea_ma_LT_soo(x,sim)
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

% arrival on mars
MJDF1_dim = MJDF1*sim.TU/(3600*24);
[kep_MA,ksun] = uplanet(MJDF1_dim, 4);
[r2, v2] = sv_from_coe(kep_MA,ksun);
% adimensionalise
r2 = r2/sim.DU;
v2 = v2/sim.DU*sim.TU;

% N REV
N_rev= round(x(3));

% vinf
vinf_mag = x(4);
alpha = x(5);
beta = x(6);

%vdep = v1 + vinf_mag*[cos(alpha)*cos(beta); sin(alpha)*cos(beta); sin(beta)] ; 
vdep = v1;
varr = v2; %since we want to rendez-vous

[output] = CW_LowLambert( r1 , r2 , vdep , varr , N_rev , TOF1 , sim.M , sim.hp , sim.kp , sim.PS , sim );

if ~isnan(output.u) % if is not nan
    T = sqrt(output.u(:,1).^2 + output.u(:,3).^2);
    CHECK_TERM_T = 0;
    if abs(max(T)) > 0.22 %0.22
        CHECK_TERM_T = 100;
    end
    obj_fun = (output.m(1) - output.m(end))/output.m(1) + CHECK_TERM_T;
else
    obj_fun = 1e+2;
end
   

end

