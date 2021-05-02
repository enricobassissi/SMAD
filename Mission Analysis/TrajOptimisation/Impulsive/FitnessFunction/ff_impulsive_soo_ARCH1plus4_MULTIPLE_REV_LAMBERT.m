function obj_fun = ff_impulsive_soo_ARCH1plus4_MULTIPLE_REV_LAMBERT(x, data, sim)

% Penalty Factor
c=1e-1;

% setting the input times
MJD01 = x(1);
TOF1 = x(2);
MJDP1 = MJD01 + TOF1; % mjd2000 passage on ast 1
TOF2 = x(3);
MJDP2 = MJDP1 + TOF2;
TOF3 = x(4);
MJDP3 = MJDP2 + TOF3; 
TOF4 =  x(5);
MJDP4 = MJDP3 + TOF4; 

% chosing which asteroid to visit
IDP = x(6); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_2 = data.PermutationMatrix(IDP,2);
asteroid_3 = data.PermutationMatrix(IDP,3);
asteroid_4 = data.PermutationMatrix(IDP,4);

% LAMBERT REVOLUTION AROUND SUN NUMBER
L_NREV_L1 = x(7);
L_NREV_L2 = x(8);
L_NREV_L3 = x(9);
L_NREV_L4 = x(10);

% the solution in case of multi revolutions is undetermined HAVE to select 
% between semimajor axis small-a option (0) or big-a option (1)
Ncase = 0; 

% Computing position and velocity of the planets in that days
% Departure from Earth
[kep_EA,ksun] = uplanet(MJD01, 3);
[rEA, vEA] = sv_from_coe(kep_EA,ksun); % km, km/s
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
ToF_EAast1_sec = TOF1*60*60*24;
ToF_ast12_sec = TOF2*60*60*24;
ToF_ast23_sec = TOF3*60*60*24;
ToF_ast34_sec = TOF4*60*60*24;

% DV calculation with lambert
% Earth -> asteroid 1
[~,~,~,~,VI_EAast1,VF_EAast1,~,~] = lambertMR(rEA,r1,ToF_EAast1_sec,ksun,0,L_NREV_L1,Ncase,0);
dv1_EAast1 = sqrt((VI_EAast1(1)-vEA(1))^2+(VI_EAast1(2)-vEA(2))^2+(VI_EAast1(3)-vEA(3))^2);
if dv1_EAast1 < sqrt(sim.C3_max) % vinf that the launcher can give max 
    dv_extra_launch = 0;
else
    % actually you would pay the difference, so we put a very big number to
    % let the optimizer to not consider this solution because too expansive
    % dv1_EAast1 - sqrt(sim.C3_max); 
    dv_extra_launch = (dv1_EAast1 - sqrt(sim.C3_max))^2; % very high number, arbitrary
end
dv2_EAast1 = sqrt((VF_EAast1(1)-v1(1))^2+(VF_EAast1(2)-v1(2))^2+(VF_EAast1(3)-v1(3))^2);

% asteroid 1 -> 2
[~,~,~,~,VI_ast12,VF_ast12,~,~] = lambertMR(r1,r2,ToF_ast12_sec,ksun,0,L_NREV_L2,Ncase,0);
dv1_ast12 = sqrt((VI_ast12(1)-v1(1))^2+(VI_ast12(2)-v1(2))^2+(VI_ast12(3)-v1(3))^2);
dv2_ast12 = sqrt((VF_ast12(1)-v2(1))^2+(VF_ast12(2)-v2(2))^2+(VF_ast12(3)-v2(3))^2);

% dV of flyby passage on asteroid 1
dv_passage_ast1 = sqrt((VI_ast12(1)-VF_EAast1(1))^2+(VI_ast12(2)-VF_EAast1(2))^2+(VI_ast12(3)-VF_EAast1(3))^2);

% kind of penalty method
% if the TOF is too little for the lambert to be solved, the algorithm put
% as result zero and this is not good
% so put now a lot on that solution to exclude it
if norm(VI_EAast1) == 0
    dv_passage_ast1 = 30;
end

% asteroid 2 -> 3
[~,~,~,~,VI_ast23,VF_ast23,~,~] = lambertMR(r2,r3,ToF_ast23_sec,ksun,0,L_NREV_L3,Ncase,0);
dv1_ast23 = sqrt((VI_ast23(1)-v2(1))^2+(VI_ast23(2)-v2(2))^2+(VI_ast23(3)-v2(3))^2);
dv2_ast23 = sqrt((VF_ast23(1)-v3(1))^2+(VF_ast23(2)-v3(2))^2+(VF_ast23(3)-v3(3))^2);

% dV of flyby passage on asteroid 2
dv_passage_ast2 = sqrt((VI_ast23(1)-VF_ast12(1))^2+(VI_ast23(2)-VF_ast12(2))^2+(VI_ast23(3)-VF_ast12(3))^2);
if norm(VI_ast23) == 0
    dv_passage_ast2 = 30;
end

% asteroid 3 -> 4
[~,~,~,~,VI_ast34,VF_ast34,~,~] = lambertMR(r3,r4,ToF_ast34_sec,ksun,0,L_NREV_L4,Ncase,0);
dv1_ast34 = sqrt((VI_ast34(1)-v3(1))^2+(VI_ast34(2)-v3(2))^2+(VI_ast34(3)-v3(3))^2);
dv2_ast34 = sqrt((VF_ast34(1)-v4(1))^2+(VF_ast34(2)-v4(2))^2+(VF_ast34(3)-v4(3))^2);

% dV of flyby passage on asteroid 3
dv_passage_ast3 = sqrt((VI_ast34(1)-VF_ast23(1))^2+(VI_ast34(2)-VF_ast23(2))^2+(VI_ast34(3)-VF_ast23(3))^2);
if norm(VI_ast34) == 0
    dv_passage_ast3 = 30;
end

% if the last dv is a flyby it can go wherever it wants
% obj_fun = dv_extra_launch + dv_passage_ast1 + dv_passage_ast2 + dv_passage_ast3; 

% penalty factor on the rel vel?
obj_fun = dv_extra_launch + dv_passage_ast1 + dv_passage_ast2 + dv_passage_ast3 + ...
    c*dv2_EAast1 + c*dv2_ast12 + c*dv2_ast23 + c*dv2_ast34; 
% else if the last dv is a rendezvous with the last object
% obj_fun = dv_extra_launch + dv_passage_ast1 + dv_passage_ast2 + dv_passage_ast3 + dv2_ast34; 

% obj_fun = TOF1+TOF2+TOF3+TOF4;

end

