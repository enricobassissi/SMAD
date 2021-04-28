function [delta_vtot, time1, time2, time3, delta, hga, norm_deltaV1, norm_deltaV3, norm_deltav_ga,e1,e2,a1,a2,norm_vinf_m,norm_vinf_p,deltav_flyby, TOF]=mer_jup_nep_ga(time)

% mer_jup_nep_ga fitness function for optimizing total velocity variation
% with hybrid method of ga and fmincon

% PROTOTYPE:
% [delta_vtot, time1, time2, time3] = mer_jup_nep_ga(time)

%INPUT:
% time        [1x3]         Normalized time variables           [-]

% OUTPUT:
% delta_vtot  [1x1]         Total velocity variation modulus    [Km/sec]
% time1       [1x1]         Time of departure in MJD2000        [days]
% time2       [1x1]         Time of flyby in MJD2000            [days]
% time3       [1x1]         Time of arrival in MJD2000          [days]

% CONTRIBUTORS:
% Berto Marina
% De Luca Maria Alessandra
% Lopes Luiz Felipe
% Sirani Samuele

% VERSIONS
%  2019-12-03 First version
%  2019-12-28 Second version
%  2020-01-18 Third version


%% Acceptable range of variables

date_a = [2028 1 1 0 0 0];
date_b = [2029 1 1 0 0 0];
mjd2000_a = date2mjd2000(date_a); 
mjd2000_b=date2mjd2000(date_b);

%flyby
date_c = [2029 1 1 0 0 0]; 
date_d = [2030 1 1 0 0 0]; 
mjd2000_c = date2mjd2000(date_c);
mjd2000_d = date2mjd2000(date_d);

%arrival
date_e = [2059 1 1 0 0 0];
date_f = [2060 1 1 0 0 0];
mjd2000_e = date2mjd2000(date_e);
mjd2000_f = date2mjd2000(date_f);

%% Normalisation of variables

time1 = time(1)*(mjd2000_b-mjd2000_a)+mjd2000_a; %time departure
time2= time(2)*(mjd2000_d-mjd2000_c)+mjd2000_c; %time flyby
time3= time(3)*(mjd2000_f-mjd2000_e)+mjd2000_e; %time arrival

%% Calculation

%kep:[a e i Om om theta]
[kep1,ksun] = uplanet(time1,1); %Mercury(departure)
[rvett1, vvett1] = kep2car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),ksun);
[kep2,ksun] = uplanet(time2,5); %Jupiter(arrival and departure)
[rvett2, vvett2] = kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),ksun);
TOF1=(time2-time1)*24*3600;
if TOF1<0
   TOF1=NaN;
end
%%%%%%%%%%%%%% first leg
MU = ksun;
[~,~,~,~,VI1,VF1,~,~] = lambertMR(rvett1,rvett2,TOF1,MU,0,0,0,0);
deltaV1=sqrt((VI1(1)-vvett1(1))^2+(VI1(2)-vvett1(2))^2+(VI1(3)-vvett1(3))^2);
if deltaV1>25
   deltaV1=NaN;
end

[kep3,ksun] = uplanet(time3,8); %Neptune(arrival)
[rvett3, vvett3] = kep2car(kep3(1),kep3(2),kep3(3),kep3(4),kep3(5),kep3(6),ksun);
TOF2=(time3-time2)*24*3600;
if TOF2<0
   TOF2=NaN;
end
%%%%%%%%%% second leg
[~,~,~,~,VI2,VF2,~,~] = lambertMR(rvett2,rvett3,TOF2,MU,0,0,0,0);
deltaV3=sqrt((VF2(1)-vvett3(1))^2+(VF2(2)-vvett3(2))^2+(VF2(3)-vvett3(3))^2);
if deltaV3>25
    deltaV3=NaN;
end
deltav_flyby=norm(VI2-VF1);

%%%%%%%%%% flyby

if time2 > time1 &&  time3 > time2
    V_m=VF1; %V_minus
    V_p=VI2; %V_plus
    mu_J= astroConstants(15);
    vinf_p=V_p-vvett2';
    vinf_m=V_m-vvett2';
    dvinf=vinf_p-vinf_m;
    norm_dvinf=norm(dvinf);
    delta=acos((-(norm_dvinf)^2+norm(vinf_p)^2+norm(vinf_m)^2)/(2*norm(vinf_p)*norm(vinf_m))); %Carnot theorem
    em =@(x) 1+ (x * norm(vinf_m)^2)/mu_J;
    ep  = @(x) 1 + (x * norm(vinf_p)^2)/mu_J;
    fun = @(x) delta - asin(1/em(x)) - asin(1/ep(x));
    rp_star=astroConstants(25); %Jupiter radius
    options=optimoptions('fsolve', 'TolFun', 1e-13, 'TolX', 1e-13,'Display','off');
    rp=fsolve(fun,rp_star, options); %distance from centre of Jupiter
    e1=1+(rp*dot(vinf_m,vinf_m)/mu_J);
    e2=1+(rp*dot(vinf_p,vinf_p)/mu_J);
    a1=rp/(1-e1);
    a2=rp/(1-e2);
    hga=rp-rp_star;
    vpm = [0;sqrt(2*(norm(vinf_m)^2/2 + mu_J/norm(rp)));0];
    vpp = [0;sqrt(2*(norm(vinf_p)^2/2 + mu_J/norm(rp)));0];
    deltav_ga = sqrt((vpm(1)-vpp(1))^2 +(vpm(2)-vpp(2))^2 + (vpm(3)-vpp(3))^2);      
    % Compute Deltavtot mission
    delta_vtot = abs(deltaV1) + abs(deltav_ga) + abs(deltaV3);
    if  hga > 0 
        delta_vtot = abs(deltaV1) + abs(deltav_ga) + abs(deltaV3);
    else
        delta_vtot = NaN;
    end               
else
    delta_vtot = NaN;
    deltav_ga = NaN;
    deltav2 = NaN;
    TOF2=NaN;
end

norm_deltaV1=norm(deltaV1);
norm_deltaV3=norm(deltaV3);
norm_deltav_ga=norm(deltav_ga);
norm_vinf_m=norm(vinf_m);
norm_vinf_p=norm(vinf_p);
TOF=TOF1+TOF2;

return