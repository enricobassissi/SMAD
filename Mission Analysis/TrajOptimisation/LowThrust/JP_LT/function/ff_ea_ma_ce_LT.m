function obj_fun = ff_ea_ma_ce_LT(x,sim)

% optimization vector
MJD01 = x(1);
TOF1  = x(2);
N_rev = round(x(3));
TOF2  = x(4);
N_rev2 = round(x(5));
rp = x(6);

% times
MJD02 = MJD01 + TOF1;
MJD03 = MJD02 + TOF2;

% Computing position and velocity of the planets in that days
% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep1,ksun] = uplanet(MJD01_dim, sim.body(1));
[r1, v1] = sv_from_coe(kep1,ksun);
% adimensionalise
r1 = r1/sim.DU;
v1 = v1/sim.DU*sim.TU;

% Flyby on Mars
mu2 = astroConstants(14);

MJD02_dim = MJD02*sim.TU/(3600*24);
[kep2,ksun] = uplanet(MJD02_dim, sim.body(2));
[r2, v2] = sv_from_coe(kep2,ksun);
% adimensionalise
r2 = r2/sim.DU;
v2 = v2/sim.DU*sim.TU;

% Arrival on Ceres
MJD03_dim = MJD03*sim.TU/(3600*24);
[kep3,ksun] = uplanet(MJD03_dim, sim.body(3));
[r3, v3] = sv_from_coe(kep3,ksun);
% adimensionalise
r3 = r3/sim.DU;
v3 = v3/sim.DU*sim.TU;

% Vasile paper
beta = 2*asin*(mu2/(v2^2*rp + mu2));

% vinfinito all'arrivo a marte?
% ruoto vinfinito dell'angolo beta

% constraints su rp 
% Delta-V slide prof 21 x

[output2] = CW_LowLambert( r2 , r3 , v2 , v3 , N_rev2 , TOF2 , sim.M , sim.hp , sim.kp , sim.PS , sim ); 

[output]  = CW_LowLambert( r1 , r2 , v1 , v2 , N_rev , TOF1 , output2.m, sim.hp , sim.kp , sim.PS , sim );

% Heliocentric velocity of the s/c arriving to body 2
vs_1 = output.v_r.*cos(output.w);
vs_2 = output.v_r.*sin(output.w);
vs_3 = output.v_z;

v_inf = [vs_1; vs_2; vs_3] - v2;
v_inf_abs = abs(v_inf);

% delta
delta = asin(mu2/(mu2 + rp*v_inf_abs^2));


obj_fun(1) = TOF1;

obj_fun(2) = TOF2;
   
obj_fun(3) = (output.m(1) - output.m(end))/output.m(1); % la massa è dimensionale

end

