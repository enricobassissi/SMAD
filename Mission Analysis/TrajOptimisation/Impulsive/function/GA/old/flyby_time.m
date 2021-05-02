function [tfb,t_m,t_p]=flyby_time(planet,a_m,e_m,a_p,e_p)
%{
CONTRIBUTORS: Apparenza Lucia
              Bassissi Enrico
              Di Trocchio Marco
              Lane John
%}
M_SUN=1.9885e30; %kg

switch planet
    case 'mars'
        rp_MA=2.067e8; %km
        ra_MA=2.492e8; %km
        a_pl=(rp_MA+ra_MA)/2; %km
        mu_pl=astroConstants(14); %[km^3/s^2]
        m_pl=6.4171e23; %kg
    case 'earth'
        rp_E=1.47095e8; %km
        ra_E=1.521e8; %km
        a_pl=(rp_E+ra_E)/2; %km
        mu_pl=astroConstants(13); %[km^3/s^2]
        m_pl=5.972e24; %kg
        
end
        
% MINUS BRANCH
% semilatus rectum
p_m=a_m*(1-e_m^2);

% sphere of influence of mars
% a_MA semi major axis mars sun
% m_MA mass of mars
% M_SUN mass of the sun
% rSOI= a_pl*(m_pl/M_SUN)^(2/5); %km
%adapted rsoi
rSOI= 0.9431*a_pl*(m_pl/M_SUN)^(2/5); %km
h_p = sqrt(mu_pl * a_p * (1 - e_p^2));

theta_p = acos((((h_p^2)/(rSOI*mu_pl))-1)/e_p);

F_p = 2 * atanh(sqrt((e_p-1)/(e_p+1))*tan(theta_p/2));

M_p = e_p * sinh(F_p) - F_p;

t_p = (M_p * h_p^3) / (mu_pl^2 * (e_p^2 - 1)^1.5);


h_m = sqrt(mu_pl* a_m * (1 - e_m^2));

theta_m = acos((((h_m^2)/(rSOI*mu_pl))-1)/e_m);

F_m = 2 * atanh(sqrt((e_m-1)/(e_m+1))*tan(theta_m/2));

M_m = e_m * sinh(F_m) - F_m;

t_m = (M_m * h_m^3) / (mu_pl^2 * (e_m^2 - 1)^1.5);

tfb = t_m + t_p; %in seconds

% tfb_d = (t_m + t_p)/(24*60*60); %in days

end
 

 