R = MA.SC1.uniform.R; % km
time = MA.SC1.uniform.time;
time_mjd2000 = MA.SC1.mjd2000_dep + time./86400;
for i = 1:length(time)
    [kep_EA, muSun] = uplanet(time_mjd2000(i),3);
    [r_EA(i,:), v_EA(i,:)] = sv_from_coe(kep_EA, muSun); % km, km/s
end

R_EA_norm = vecnorm(r_EA,2,2); % km

Dist_Rel = vecnorm(R - r_EA,2,2);

figure()
plot(time_mjd2000,Dist_Rel);
xlabel('mjd2000'); ylabel('Dist SC - Earth [km]');