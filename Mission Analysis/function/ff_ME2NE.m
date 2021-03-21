function dvtot = ff_ME2NE(T)
%{
UNTITLED4 Summary of this function goes here
t1 Departure date in mJ2000
t2 Arrival date in mJ2000

CONTRIBUTORS: Apparenza Lucia
              Bassissi Enrico
              Di Trocchio Marco
              Lane John
%}
% Computing position and velocity of the planets in that days
[kep_ME,ksun] = uplanet(T(1), 1);
[r1, v1] = sv_from_coe(kep_ME,ksun);

[kep_MA,ksun] = uplanet(T(2), 4);
[r2, v2] = sv_from_coe(kep_MA,ksun);

[kep_NE,ksun] = uplanet(T(3), 8);
[r3, v3] = sv_from_coe(kep_NE,ksun); 

R_MA=astroConstants(24);  % Mars radius 
h_lim=200;
shift=200;
mu_MA=astroConstants(14); % Gravitational constant of mars

% Converting mJ2000 in seconds
t1_sec=T(1).*60.*60.*24;
t2_sec=T(2).*60.*60.*24;
t3_sec=T(3).*60.*60.*24;

% DV calculation with lambert
[dvtot12,~,VF12]=lambert_solver_flyby( r1, r2, v1, v2, t1_sec, t2_sec, ksun, 1e2);
[dvtot23,VI23,~]=lambert_solver_flyby2( r2, r3, VF12, v3, t2_sec, t3_sec, ksun, 1e2);
%[dvGA,rp] = PGA_solver (VF12,VI23,v2,mu_MA,h_lim+shift,R_MA);
%[dvGA,rp] = poweredGA (VF12,VI23,R_MA,v2,mu_MA);
%[~,dvGA,~]=gravity_assist(VF12,VI23,v2,mu_MA,R_MA,h_lim);
% [dvGA,~]=PGA(VF12,VI23,v2,R_MA,mu_MA,h_lim);
[dvGA,~]=PGA_ga(VF12,VI23,v2,R_MA,mu_MA,h_lim);
% dvtot(1)=dvtot12;
% dvtot(2)=dvtot23;
% dvtot(3)=dvGA;
dvtot=dvtot12+dvtot23+dvGA;

end

