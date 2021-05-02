function obj_fun = ff_impulsive_soo_ARCH1plus4_mass(x, data, sim)
% setting the input times
MJD01 = x(1);
TOF0 = x(2);
MJDGA = MJD01 + TOF0;
TOF1 = x(3);
MJDP1 = MJDGA + TOF1; % mjd2000 passage on ast 1
TOF2 = x(4);
MJDP2 = MJDP1 + TOF2;
TOF3 = x(5);
MJDP3 = MJDP2 + TOF3; 
TOF4 =  x(6);
MJDP4 = MJDP3 + TOF4; 

% chosing which asteroid to visit
IDP = round(x(7)); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_2 = data.PermutationMatrix(IDP,2);
asteroid_3 = data.PermutationMatrix(IDP,3);
asteroid_4 = data.PermutationMatrix(IDP,4);

% Computing position and velocity of the planets in that days
% Departure from Earth
[kep_EA,ksun] = uplanet(MJD01, 3);
[rEA, vEA] = sv_from_coe(kep_EA,ksun); % km, km/s
% Earth GA
[kep_GA,~] = uplanet(MJDGA, 3);
[rGA, vGA] = sv_from_coe(kep_GA,ksun); % km, km/s
% passage at 1st ast
[kep_ast_1] = uNEO2(MJDP1,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1, v1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
% passage at 2nd asteroid 
[kep_ast_2] = uNEO2(MJDP2,asteroid_2,data);
[r2, v2] = sv_from_coe(kep_ast_2,ksun); % km, km/s
% passage at 3rd asteroid 
[kep_ast_3] = uNEO2(MJDP3,asteroid_3,data);
[r3, v3] = sv_from_coe(kep_ast_3,ksun); % km, km/s
% arrival at 4th asteroid 
[kep_ast_4] = uNEO2(MJDP4,asteroid_4,data);
[r4, v4] = sv_from_coe(kep_ast_4,ksun); % km, km/s

% Converting TOFs in seconds for Lambert
ToF_EAGA_sec = TOF0*60*60*24;
ToF_GAast1_sec = TOF1*60*60*24;
ToF_ast12_sec = TOF2*60*60*24;
ToF_ast23_sec = TOF3*60*60*24;
ToF_ast34_sec = TOF4*60*60*24;

% DV calculation with lambert
% Earth -> asteroid 1
% with Nrev = 1 and Ncase = 0, it find only sol tha follow the earth
[~,~,~,~,VI_EAGA,VF_EAGA,~,~] = lambertMR(rEA,rGA,ToF_EAGA_sec,ksun,0,0,0,0);

dv1_EAGA = sqrt((VI_EAGA(1)-vEA(1))^2+(VI_EAGA(2)-vEA(2))^2+(VI_EAGA(3)-vEA(3))^2);
if dv1_EAGA < sqrt(sim.C3_max) % vinf that the launcher can give max 
    dv_extra_launch = 0;
else
    % actually you would pay the difference, so we put a very big number to
    % let the optimizer to not consider this solution because too expansive
    % dv_extra_launch = dv1_EAast1 - sqrt(sim.C3_max); 
    % dv_extra_launch = 20; % very high number, arbitrary
    c1 = 5; % penalty factor for dv_extra_launch
    dv_extra_launch = c1*(dv1_EAGA - sqrt(sim.C3_max))^2; % penalty like, but not discard a priori
end
% dv2_EAGA = sqrt((VF_EAGA(1)-vGA(1))^2+(VF_EAGA(2)-vGA(2))^2+(VF_EAGA(3)-vGA(3))^2);

[~,~,~,~,VI_GAast1,VF_GAast1,~,~] = lambertMR(rGA,r1,ToF_GAast1_sec,ksun,0,0,0,0);
dv2_GAast1 = sqrt((VF_GAast1(1)-v1(1))^2+(VF_GAast1(2)-v1(2))^2+(VF_GAast1(3)-v1(3))^2);

% astroConstants(23) = Radius_Earth, km
% astroConstants(13) = muEarth, km^3/s^2
delta_V_p = flyby(astroConstants(23), astroConstants(13), MJDGA, VF_EAGA, VI_GAast1);

if strcmp(string(delta_V_p), 'Not found')
    delta_V_p = 30;
end

% asteroid 1 -> 2
[~,~,~,~,VI_ast12,VF_ast12,~,~] = lambertMR(r1,r2,ToF_ast12_sec,ksun,0,0,0,0);
% dv1_ast12 = sqrt((VI_ast12(1)-v1(1))^2+(VI_ast12(2)-v1(2))^2+(VI_ast12(3)-v1(3))^2);
dv2_ast12 = sqrt((VF_ast12(1)-v2(1))^2+(VF_ast12(2)-v2(2))^2+(VF_ast12(3)-v2(3))^2);

% dV of flyby passage on asteroid 1
dv_passage_ast1 = sqrt((VI_ast12(1)-VF_GAast1(1))^2+(VI_ast12(2)-VF_GAast1(2))^2+(VI_ast12(3)-VF_GAast1(3))^2);
    
% asteroid 2 -> 3
[~,~,~,~,VI_ast23,VF_ast23,~,~] = lambertMR(r2,r3,ToF_ast23_sec,ksun,0,0,0,0);
% dv1_ast23 = sqrt((VI_ast23(1)-v2(1))^2+(VI_ast23(2)-v2(2))^2+(VI_ast23(3)-v2(3))^2);
dv2_ast23 = sqrt((VF_ast23(1)-v3(1))^2+(VF_ast23(2)-v3(2))^2+(VF_ast23(3)-v3(3))^2);

% dV of flyby passage on asteroid 2
dv_passage_ast2 = sqrt((VI_ast23(1)-VF_ast12(1))^2+(VI_ast23(2)-VF_ast12(2))^2+(VI_ast23(3)-VF_ast12(3))^2);

% asteroid 3 -> 4
[~,~,~,~,VI_ast34,VF_ast34,~,~] = lambertMR(r3,r4,ToF_ast34_sec,ksun,0,0,0,0);
% dv1_ast34 = sqrt((VI_ast34(1)-v3(1))^2+(VI_ast34(2)-v3(2))^2+(VI_ast34(3)-v3(3))^2);
dv2_ast34 = sqrt((VF_ast34(1)-v4(1))^2+(VF_ast34(2)-v4(2))^2+(VF_ast34(3)-v4(3))^2);

% dV of flyby passage on asteroid 3
dv_passage_ast3 = sqrt((VI_ast34(1)-VF_ast23(1))^2+(VI_ast34(2)-VF_ast23(2))^2+(VI_ast34(3)-VF_ast23(3))^2);

% if the last dv is a flyby it can go wherever it wants
% obj_fun = dv_extra_launch + dv_passage_ast1 + dv_passage_ast2 + dv_passage_ast3 + delta_V_p; 

% else if the last dv is a rendezvous with the last object
% obj_fun = dv_extra_launch + dv_passage_ast1 + dv_passage_ast2 + dv_passage_ast3 + dv2_ast34; 

% % Penalty Factor
% c=1e-1;
% % penalty factor on the rel vel?
% obj_fun = dv_extra_launch + dv_passage_ast1 + dv_passage_ast2 + dv_passage_ast3 + ...
%     c*dv2_GAast1 + c*dv2_ast12 + c*dv2_ast23 + c*dv2_ast34; 

% Mass Calculations
g0 = sim.g0; %km/s^2
Isp_mother=sim.Isp_mother; %s
Isp_lander=sim.Isp_lander; %s
dry_mass=sim.dry_mass; %kg
dry_mass_lander=sim.dry_mass_lander; %kg
% m_wet/m_dry = exp(dv/(Isp*g0))
% mass of the lander required to pay to cancel dVrel at each flyby passage
% dry mass lander and the propellant related to that dV
mass_lander=[dry_mass_lander*exp(dv2_GAast1/(Isp_lander*g0)), dry_mass_lander*exp(dv2_ast12/(Isp_lander*g0)),...
             dry_mass_lander*exp(dv2_ast23/(Isp_lander*g0)), dry_mass_lander*exp(dv2_ast34/(Isp_lander*g0))];
% Back interpolate the mass of the overall spacecraft from end of mission,
% with an end mass of sim.dry_mass and we build back the initial wet mass
% at each asteroid encounter, the mothercraft expell a lander
bef_ast4=dry_mass+mass_lander(4);
% at each flyby the sc expell mass to perform the impulsive dV
aft_ast3=bef_ast4*exp((dv_passage_ast3)/(Isp_mother*g0));
bef_ast3=aft_ast3+mass_lander(3);
aft_ast2=bef_ast3*exp((dv_passage_ast2)/(Isp_mother*g0));
bef_ast2=aft_ast2+mass_lander(2);
aft_ast1=bef_ast2*exp((dv_passage_ast1)/(Isp_mother*g0));
bef_ast1=aft_ast1+mass_lander(1);
% mass expelled from powered gravity assist on the earth
aft_ga = bef_ast1*exp((delta_V_p)/(Isp_mother*g0));
% mass used due to the assisted insertion in the 1st interplanetary leg,
% other then the sqrt(C3) of the launcher, if any
wet_mass=aft_ga*exp((dv_extra_launch)/(Isp_mother*g0));

obj_fun = wet_mass;

end

