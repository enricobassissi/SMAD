function [ a ] = find_a_prof(x,psi,TOF,sim,s1,s1_x,s1_xx,s2,s2_x,s2_xx,delta1,delta1_x,delta1_xx,delta2,delta2_x,delta2_xx)

n_sol = sim.n_sol;

% Ds Ds_x  Ds_xx 
Ds     = s2 - s1;
Ds_x   = s2_x - s1_x;
Ds_xx  = s2_xx - s1_xx;


% Ddelta Ddelta_x Ddelta_xx Ddelta_xxx 
Ddelta     = delta2 - delta1 ;
Ddelta_x   = delta2_x - delta1_x;
Ddelta_xx  = delta2_xx - delta1_xx;

a = 0;
% xi = xi(x,a)
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
x_t_2 = sim.mu*r./(s.^3.*(r*psi^2 - r_xx + 2*r_x.^2./r));  % \dot{x}^2
x_t   = sqrt(x_t_2); % \dot{x}

% step zero di newton
% dx = x(2) - x(1);
% %dx = 1e-8;
% j = 3:2:sim.n_sol-2;
% k = 2:2:sim.n_sol-1;
% 
% dTOF1 = 1./x_t(1);
% dTOFend = 1./x_t(sim.n_sol);
% 
% dTOFk = 1./x_t(k); %dTOFm
% dTOFj = 1./x_t(j);
% 
% TOFc1 = (dTOF1 + dTOFend)/6 + sum(dTOFj)/3;
% TOFc1 = 2*dx*(TOFc1 + 2/3*sum(dTOFk))
integrand = 1./x_t;
TOFc1 = trapz(x,integrand);

error = (TOF - TOFc1)/TOF;

a_old = a;
a = 0.00001;
count = 0;


while abs(error) > 1e-6 && count <=100
    
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
    r_xx  = - (2*s_x.*delta_x + s.*delta_xx).*sin(delta) + (s_xx - s.*delta_x.^2).*cos(delta);


    % Derivative of x
    x_t_2 = sim.mu*r./(s.^3.*(r*psi^2 - r_xx + 2*r_x.^2./r)); 
    x_t   = sqrt(x_t_2);
    
%     dTOF1 = 1./x_t(1); % first term 1/\dot{x}
%     dTOFend = 1./x_t(sim.n_sol); % last term
%     
%     dTOFk = 1./x_t(k); %dTOFm
%     dTOFj = 1./x_t(j);
%     
%     TOFc = (dTOF1 + dTOFend)/6 + sum(dTOFj)/3;
%     TOFc = 2*dx*(TOFc + 2/3*sum(dTOFk))
    integrand = 1./x_t;
    TOFc = trapz(x,integrand);
    error_old = error;
    error = (TOF - TOFc)/TOF;
   % d_error = (error - error_old)/(a - a_old)
    d_error = (error - error_old)/dx; %%% giusto????
    
    % newton
    a_old = a;
    a = a_old - error/d_error; 
    
    count = count + 1
end
find_d
% if count >= 99
%     T = 9999999 + 0*x;

end