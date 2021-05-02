xm   = 0.5*(x(1:end-1) + x(2:end));

psim = psi; %%% controlla

% Interpolating function 
xi     = @(a) a*x^8 - (20 + 4*a)*x^7 + (70 + 6*a)*x^6 - (84 + 4*a)*x^5 + (35 + a)*x^4;
xi_x   = @(a) 8*a*x^7 - 7*(20+4*a)*x^6 + 6*(70 + 6*a)*x^5 - 5*(84 + 4*a)*x^4 + 4*(35 + a)*x^3;
xi_xx  = @(a) 56*a*x^6 - 42*(20+4*a)*x^5 + 30*(70 + 6*a)*x^4 - 20*(84 + 4*a)*x^3 + 12*(35 + a)*x^2;

xim     = @(a) a*xm^8 - (20 + 4*a)*xm^7 + (70 + 6*a)*xm^6 - (84 + 4*a)*xm^5 + (35 + a)*xm^4;
xim_x   = @(a) 8*a*xm^7 - 7*(20+4*a)*xm^6 + 6*(70 + 6*a)*xm^5 - 5*(84 + 4*a)*xm^4 + 4*(35 + a)*xm^3;
xim_xx  = @(a) 56*a*xm^6 - 42*(20+4*a)*xm^5 + 30*(70 + 6*a)*xm^4 - 20*(84 + 4*a)*xm^3 + 12*(35 + a)*xm^2;

% Ds Ds_x  Ds_xx %%%% controlla
Ds     = s2 - s1;
Ds_x   = s2_x - s1_x;
Ds_xx  = s2_xx - s1_xx;

% Ddelta Ddelta_x Ddelta_xx Ddelta_xxx %%%% controlla
Ddelta     = delta2 - delta1 ;
Ddelta_x   = delta2_x - delta1_x;
Ddelta_xx  = delta2_xx - delta1_xx;

% Sun distance s(x) and declination delta(x) + derivatives
s     = @(a) (s2 - s1) * xi(a) + s1;
delta = @(a) (delta2 - delta1) * xi(a) + delta1; 

s_x   = @(a) Ds_x*xi(a) + xi_x(a) *Ds + s1_x;
s_xx  = @(a) Ds_xx*xi(a) + xi_xx(a)*Ds + 2*Ds_x*xi_x(a) + s1_xx;

delta_x   = @(a) Ddelta_x*xi(a) + xi_x(a)*Ddelta + delta1_x;
delta_xx  = @(a) Ddelta_xx*xi(a) + xi_xx(a)*Ddelta + 2*Ddelta_x*xi_x(a) + delta1_xx;

sm     = @(a) (s2 - s1) * xim(a) + s1;
deltam = @(a) (delta2 - delta1) * xim(a) + delta1; 

sm_x   = @(a) Ds_x*xim(a) + xim_x(a) *Ds + s1_x;
sm_xx  = @(a) Ds_xx*xim(a) + xim_xx(a)*Ds + 2*Ds_x*xim_x(a) + s1_xx;

deltam_x   = @(a) Ddelta_x*xim(a) + xim_x(a)*Ddelta + delta1_x;
deltam_xx  = @(a) Ddelta_xx*xim(a) + xim_xx(a)*Ddelta + 2*Ddelta_x*xim_x(a) + delta1_xx;

% r and its derivatives 
r     = @(a) s(a)*cos(delta(a));
r_x   = @(a) - s(a)*delta_x(a) + s_x(a) * cos(delta(a));
r_xx  = @(a) - (2*s_x(a)*delta_x(a) + s_x(a)*delta_xx(a))*sin(delta(a)) + (s_xx(a) - s(a)*delta_x(a)^2)*cos(delta(a));

rm     = @(a) sm(a)*cos(deltam(a));
rm_x   = @(a) - sm(a)*deltam_x(a) + sm_x(a) * cos(deltam(a));
rm_xx  = @(a) - (2*sm_x(a)*deltam_x(a) + sm_x(a)*deltam_xx(a))*sin(deltam(a)) + (sm_xx(a) - sm(a)*deltam_x(a)^2)*cos(deltam(a));

% Derivative of x
x_t_2 = @(a) sim.mu*r(a)/(s(a)^3*(r(a)*psi^2 - r_xx(a) + 2*r_x(a)^2/r(a)));
x_t   = @(a) sqrt(x_t_2(a)); 

xm_t_2 = @(a) sim.mu*rm(a)/(sm(a)^3*(rm(a)*psim^2 - rm_xx(a) + 2*rm_x(a)^2/rm(a)));
xm_t   = @(a) sqrt(xm_t_2(a)); 

% dTOF at each theta 
dTOF  = 1./x_t(a);
dTOFm = 1./xm_t(a);

% step size
h = x(2)- x(1);

% Cavalieri-Simpson method
I = 0;
for i = 2:n_sol
    I = I + dTOF(i-1) + 4*dTOFm(i-1) + dTOF(i);
end

I = I*h/6;

% Computation of TOF error (residual)
error_TOF =TOF-I;
