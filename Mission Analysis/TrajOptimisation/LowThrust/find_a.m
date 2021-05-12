function [ error_TOF ] = find_a(a,x,psi,TOF,sim,s1,s1_x,s1_xx,s2,s2_x,s2_xx,delta1,delta1_x,delta1_xx,delta2,delta2_x,delta2_xx)

n_sol = sim.n_sol;

% Ds Ds_x  Ds_xx 
Ds     = s2 - s1;
Ds_x   = s2_x - s1_x;
Ds_xx  = s2_xx - s1_xx;


% Ddelta Ddelta_x Ddelta_xx Ddelta_xxx 
Ddelta     = delta2 - delta1 ;
Ddelta_x   = delta2_x - delta1_x;
Ddelta_xx  = delta2_xx - delta1_xx;

xm   = 0.5*(x(1:end-1) + x(2:end));


xi     = a*x.^8 - (20 + 4*a)*x.^7 + (70 + 6*a)*x.^6 - (84 + 4*a)*x.^5 + (35 + a)*x.^4;
xi_x   = 8*a*x.^7 - 7*(20+4*a)*x.^6 + 6*(70 + 6*a)*x.^5 - 5*(84 + 4*a)*x.^4 + 4*(35 + a)*x.^3;
xi_xx  = 56*a*x.^6 - 42*(20+4*a)*x.^5 + 30*(70 + 6*a)*x.^4 - 20*(84 + 4*a)*x.^3 + 12*(35 + a)*x.^2;

% Sun distance s(x) and declination delta(x) + derivatives
s     =  (s2 - s1) .* xi + s1;
delta =  (delta2 - delta1) .* xi + delta1; 

s_x   =  Ds_x.*xi + xi_x .*Ds + s1_x;
s_xx  =  Ds_xx.*xi + xi_xx.*Ds + 2*Ds_x.*xi_x + s1_xx;

delta_x   = Ddelta_x.*xi + xi_x.*Ddelta + delta1_x;
delta_xx  = Ddelta_xx.*xi + xi_xx.*Ddelta + 2*Ddelta_x.*xi_x + delta1_xx;

% r and its derivatives 
r     = s.*cos(delta);
r_x   = - s.*delta_x.*sin(delta) + s_x .* cos(delta);
r_xx  = - (2*s_x.*delta_x + s_x.*delta_xx).*sin(delta) + (s_xx - s.*delta_x.^2).*cos(delta);

% Derivative of x
x_t_2 = sim.mu*r./(s.^3.*(r*psi^2 - r_xx + 2*r_x.^2./r));  
x_t   = sqrt(x_t_2);

xm_t = 0.5*(x_t(1:end-1) + x_t(2:end));

% ----------------------------------------------------------------------- %
% MIA 1
% dTOF at each theta 
dTOF  = 1./x_t;
dTOFm = 1./xm_t;

% step size
h = x(2)- x(1);
%h = 1e-8;

% Cavalieri-Simpson method
I = 0;
for i = 2:n_sol
    I = I + dTOF(i-1) + 4*dTOFm(i-1) + dTOF(i);
end

I = I*h/6;
error_TOF = (TOF-I);

% %------------------------------------------------------------------------%
% % MIA 2 :
% %dtheta = theta(2)- theta(1);
% dtheta = 1e-8;
% dtheta6 = dtheta./6;
% 
% TOF11    = theta;
% TOF11(1) = 0; % Initialization
% dTOF  = 1./x_t;
% 
% for i = 2:n_sol-1
%     TOF11(i) = TOF11(i-1) + dtheta6 * (dTOF(i-1) + 4*dTOF(i) + dTOF(i+1) ) ; %% RK4 (?)
% end
% TOF11(n_sol) = TOF11(n_sol-1) + dTOF(n_sol-1)*dtheta ;
% 
% %Computation of TOF error (residual)
% error_TOF = (TOF-TOF11); %% abs?

end