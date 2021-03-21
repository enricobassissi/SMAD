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


[kep_NE,ksun] = uplanet(T(2), 8);
[r2, v2] = sv_from_coe(kep_NE,ksun); 


% Converting mJ2000 in seconds
t1_sec=T(1).*60.*60.*24;
t2_sec=T(2).*60.*60.*24;

% DV calculation with lambert
[dvtot,~,~]=lambert_solver_flyby( r1, r2, v1, v2, t1_sec, t2_sec, ksun, 1e2);


end

