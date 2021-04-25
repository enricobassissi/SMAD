function [dvGA,r_peri,delta,tfb]=PGA(Vminus,Vplus,Vplanet,Rplanet,mu_planet,h_lim_P,planet)
%{
Function used to compute the delta v required  by the powered gravity
assist manoeuvre and the parameters of the arrival and departure
hyperbola. 

INPUT 
V_inf_m       [1x3] :Heliocentric velocity of the arrival branch [km/s]
V_inf_p       [1x3] :Heliocentric velocity of the departure branch [km/s]
ibody_fb      [1] : Integer defining the body near which the flyby is performed
date          [1] : Date of the flyby, expressed in modified Julian days

          
OUTPUT 
r_peri                radius of the pericenter [km]
e_1                   eccentricity pf the incoming branch
e_2                   eccentricity of the outgoing branch
vp_1[1]               Velocity at the pericenter of the arrival fly-by hyperbola  [km/s].
vp_2[1]               Velocity at the pericenter of the departure fly-by hyperbola  [km/s].
dv [1]                Delta v required to perform the manoeuvre [km/s]
time [1]              Flyby Time [s]

CONTRIBUTORS: Apparenza Lucia
              Bassissi Enrico
              Di Trocchio Marco
              Lane John
%}

% Excess infinite velocity with respect to the planet 
vinf_m = Vminus - Vplanet;
vinf_p = Vplus - Vplanet;
% Determine the pericenter radius of the hyperbola and the other parameters
% using fsolve 
% Fsolve options 
delta=acos(dot(vinf_p,vinf_m)/(norm(vinf_p)*norm(vinf_m)));
options = optimoptions('fsolve','Display','none', 'OptimalityTolerance', 1e-8 ); 
% Solve the function of delta - delta_m - delta+ = 0 wrt r_peri 
fun = @(r_p) acos((dot(vinf_m,vinf_p)/(norm(vinf_m)*norm(vinf_p)))) - ...
     asin(1/(1+ (r_p*norm(vinf_m)^2/mu_planet))) - asin(1/(1+ (r_p*norm(vinf_p)^2/mu_planet)));
 
[r_peri]=fsolve(fun,Rplanet,options);
rperi_lim=Rplanet+h_lim_P;
if r_peri<=rperi_lim
   % disp('error, r_{perigee} not feasable, the s/c will crash');
   dvGA=NaN;
else
% Compute the parameters of the flyby hyperbola 
 a_m = -mu_planet/norm(vinf_m)^2;
 a_p = -mu_planet/norm(vinf_p)^2;
 vp_1 = sqrt(-mu_planet*(1/a_m - 2/r_peri));
 vp_2 = sqrt(-mu_planet*(1/a_p - 2/r_peri));
 e_m =  1 + r_peri*norm(vinf_m)^2 / mu_planet;
 e_p = 1 + r_peri*norm(vinf_p)^2 / mu_planet;
 %Delta V to be provided 
 dvGA = norm(vp_2-vp_1);
 [tfb,~,~]=flyby_time(planet,a_m,e_m,a_p,e_p);
end




end