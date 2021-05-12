function obj_fun = ff_impulsive_soo_ARCH1plus4_MARS_GA(x, data, sim)
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
% Mars GA
[kep_GA,~] = uplanet(MJDGA, 4);
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
    c_launcher = 30; % penalty factor for dv_extra_launch
    dv_extra_launch = c_launcher*(dv1_EAGA - sqrt(sim.C3_max))^2; % penalty like, but not discard a priori
end
% dv2_EAGA = sqrt((VF_EAGA(1)-vGA(1))^2+(VF_EAGA(2)-vGA(2))^2+(VF_EAGA(3)-vGA(3))^2);

[~,~,~,~,VI_GAast1,VF_GAast1,~,~] = lambertMR(rGA,r1,ToF_GAast1_sec,ksun,0,0,0,0);
dv2_GAast1 = sqrt((VF_GAast1(1)-v1(1))^2+(VF_GAast1(2)-v1(2))^2+(VF_GAast1(3)-v1(3))^2);

% astroConstants(24) = Radius_Mars, km
% astroConstants(14) = muMars, km^3/s^2
R_lim_from_planet = 200; % km, for mars it' ok
delta_v_p = flyby(astroConstants(24), astroConstants(14),R_lim_from_planet, MJDGA, VF_EAGA, VI_GAast1);

if strcmp(string(delta_v_p), 'Not found')
    % Arbitrarly big barrier number, related to no flyby solution
    delta_v_p = 200;
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
% obj_fun = dv_extra_launch + dv_passage_ast1 + dv_passage_ast2 + dv_passage_ast3 + delta_v_p; 

% else if the last dv is a rendezvous with the last object
% obj_fun = dv_extra_launch + dv_passage_ast1 + dv_passage_ast2 + dv_passage_ast3 + dv2_ast34; 

% % Penalty Factor
% c=0.05;
% % penalty factor on the rel vel?
% obj_fun = dv_extra_launch + dv_passage_ast1 + dv_passage_ast2 + dv_passage_ast3 + delta_v_p + ...
%     c*dv2_GAast1 + c*dv2_ast12 + c*dv2_ast23 + c*dv2_ast34; 

% Check of feasibility
CHECK_TERM = 0; CHECK_TERM_TOF = 0;
CHECK_TERM_A = 0; CHECK_TERM_B = 0; CHECK_TERM_C = 0; CHECK_TERM_D = 0;
% Mission TOF
tot_TOF = TOF1+TOF2+TOF3+TOF4;
if tot_TOF > 12*365
    CHECK_TERM_TOF = 100;
end
% Rel Velocity at passage
if dv2_GAast1 > 7
    CHECK_TERM_A = dv2_GAast1;
end
if dv2_ast12 > 7
    CHECK_TERM_B = dv2_ast12;
end
if dv2_ast23 > 7
    CHECK_TERM_C = dv2_ast23;
end
if dv2_ast34 > 7
    CHECK_TERM_D = dv2_ast34;
end
CHECK_TERM = CHECK_TERM_A+CHECK_TERM_B+CHECK_TERM_C+CHECK_TERM_D;
% 
% % Penalty Factor
% c_dVpass_1 = 20; c_dVpass_2 = 20; c_dVpass_3 = 20; c_dVdeltavp = 20;
% c_dVrel_1 = 0.1; c_dVrel_2 = 0.08; c_dVrel_3 = 0.06; c_dVrel_4 = 0.04;
% 
% avg_dVrel = mean([dv2_GAast1, dv2_ast12, dv2_ast23, dv2_ast34]);
% 
% % obj_fun = dv_extra_launch + c_dVpass_1*dv_passage_ast1 + c_dVdeltavp*delta_v_p + ...
% %     c_dVpass_2*dv_passage_ast2 + c_dVpass_3*dv_passage_ast3 + ...
% %     c_dVrel_1*(dv2_GAast1-avg_dVrel)^2 + c_dVrel_2*(dv2_ast12-avg_dVrel)^2 +...
% %     c_dVrel_3*(dv2_ast23-avg_dVrel)^2 + c_dVrel_4*(dv2_ast34-avg_dVrel)^2 + ...
% %     CHECK_TERM + CHECK_TERM_TOF; 
% obj_fun = dv_extra_launch + dv_passage_ast1 + delta_v_p + ...
%     dv_passage_ast2 + dv_passage_ast3 + ...
%     CHECK_TERM + CHECK_TERM_TOF;
obj_fun = dv_extra_launch + delta_v_p + dv_passage_ast1 + ...
    dv_passage_ast2 + dv_passage_ast3 + ...
	CHECK_TERM + CHECK_TERM_TOF;

end

