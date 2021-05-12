function [ error_TOF ] = find_a_old(a,x,psi,TOF, sin_alpha1, sin_alpha2, p1, f1, g1, L1,p2, f2, g2, L2,sim)


x = sim.x;
xm   = 0.5*(x(1:end-1) + x(2:end));

n_sol = sim.n_sol;

A = linspace(-500,500,1000);
%------------------------------------------------------------------------ %
% DEPARTURE ORBIT
for j = 1:length(A)
    a = A(j);
%- Angle beta1(x) and derivatives needed to compute the declination
 beta1     = acos(sin_alpha1*cos(psi*x));
 beta1_x   = psi*sin_alpha1*sin(psi*x)./sin(beta1);
 beta1_xx  = (psi^2*sin_alpha1*cos(psi*x) - cos(beta1).*beta1_x.^2)./sin(beta1);

 beta1m     = acos(sin_alpha1*cos(psi*xm));
 beta1m_x   = psi*sin_alpha1*sin(psi*xm)./sin(beta1m);
 beta1m_xx  = (psi^2*sin_alpha1*cos(psi*xm) - cos(beta1m).*beta1m_x.^2)./sin(beta1m);
 
 %- Declination above the ref. plane delta1(x) - è univoc definita?? %%%%
 delta1 = asin(sin_alpha1*sin(psi*x)./sin(beta1));
 
 delta1m = asin(sin_alpha1*sin(psi*xm)./sin(beta1m));
 
 %- Derivatives of delta1(x) wrt x : delta1_x delta1_xx 
 delta1_x   = (psi *sin_alpha1 * cos(psi*x) - beta1_x.*cos(beta1).*sin(delta1))./(cos(delta1) .* sin(beta1));
 delta1_xx  = (- psi^2 *sin_alpha1 .* sin(psi*x) ...
              +  sin(beta1).*sin(delta1).*(delta1_x.^2 + beta1_x.^2))./(cos(delta1) .* sin(beta1)) ...
              - (2*delta1_x.*beta1_x.*cos(delta1).*cos(beta1) ...
              +  beta1_xx.*sin(delta1).*cos(beta1))./(cos(delta1) .* sin(beta1));
          
 delta1m_x   = (psi *sin_alpha1 * cos(psi*xm) - beta1m_x.*cos(beta1m).*sin(delta1m))./(cos(delta1m) .* sin(beta1m));
 delta1m_xx  = (- psi^2 *sin_alpha1 .* sin(psi*xm) ...
              +  sin(beta1m).*sin(delta1m).*(delta1m_x.^2 + beta1m_x.^2))./(cos(delta1m) .* sin(beta1m)) ...
              - (2*delta1m_x.*beta1m_x.*cos(delta1m).*cos(beta1m) ...
              +  beta1m_xx.*sin(delta1m).*cos(beta1m))./(cos(delta1m) .* sin(beta1m));
  
 %- Variation of true longitude DL1(x)
 sin_DL1 = 1/sin_alpha1 * sin(delta1);
 cos_DL1 = cos(psi*x) .* cos(delta1);
 DL1 = atan2(sin_DL1,cos_DL1);
 
 sin_DL1m = 1/sin_alpha1 * sin(delta1m);
 cos_DL1m = cos(psi*xm) .* cos(delta1m);
 DL1m = atan2(sin_DL1m,cos_DL1m);
 
 %- Derivatives of DL1(x) wrt x * DL1_x DL1_xx DL1_xxx
 DL1_x    = delta1_x./(sin_alpha1.*cos(psi*x));
 DL1_xx   = (delta1_xx + psi*sin_alpha1.*sin(psi*x) .*DL1_x)./(sin_alpha1.*cos(psi*x));
                   
 DL1m_x    = delta1m_x./(sin_alpha1.*cos(psi*xm));
 DL1m_xx   = (delta1m_xx + psi*sin_alpha1.*sin(psi*xm) .*DL1m_x)./(sin_alpha1.*cos(psi*xm));
 
%------------------------------------------------------------------------ %
% ARRIVAL ORBIT

%- Angle beta2(x) and derivatives needed to compute the declination
 beta2     = acos(sin_alpha2*cos(psi*(1 -x)));
 beta2_x   = - psi*sin_alpha2*sin(psi*(1-x))./sin(beta2);
 beta2_xx  = (- psi^2*sin_alpha2.*cos(psi*(1-x)) - cos(beta2).*beta2_x.^2)./sin(beta2);
 
 beta2m     = acos(sin_alpha2*cos(psi*(1 -xm)));
 beta2m_x   = - psi*sin_alpha2*sin(psi*(1-xm))./sin(beta2m);
 beta2m_xx  = (- psi^2*sin_alpha2.*cos(psi*(1-xm)) - cos(beta2m).*beta2m_x.^2)./sin(beta2m);

 
 %- Declination above the ref. plane delta2(x) - è univoc definita?? %%%%
 delta2 = asin(sin_alpha2*sin(psi*(1-x))./sin(beta2));
 
 delta2m = asin(sin_alpha2*sin(psi*(1-xm))./sin(beta2m));
 
 %- Derivatives of delta2(x) wrt x : delta2_x delta2_xx delta2_xxx
 delta2_x   = (- psi *sin_alpha2 * cos(psi*(1-x)) - beta2_x.*cos(beta2).*sin(delta2))./(cos(delta2) .* sin(beta2));
 delta2_xx  = (- psi^2 *sin_alpha2 * sin(psi*(1-x)) ...
              +  sin(beta2).*sin(delta2).*(delta2_x.^2 + beta2_x.^2))./(cos(delta2) .* sin(beta2)) ...
              - (2*delta2_x.*beta2_x .*cos(delta2) .*cos(beta2) ...
              +  beta2_xx .*sin(delta2) .*cos(beta2))./(cos(delta2) .* sin(beta2));
          
 delta2m_x   = (- psi *sin_alpha2 * cos(psi*(1-xm)) - beta2m_x.*cos(beta2m).*sin(delta2m))./(cos(delta2m) .* sin(beta2m));
 delta2m_xx  = (- psi^2 *sin_alpha2 * sin(psi*(1-xm)) ...
              +  sin(beta2m).*sin(delta2m).*(delta2m_x.^2 + beta2m_x.^2))./(cos(delta2m) .* sin(beta2m)) ...
              - (2*delta2m_x.*beta2m_x .*cos(delta2m) .*cos(beta2m) ...
              +  beta2m_xx .*sin(delta2m) .*cos(beta2m))./(cos(delta2m) .* sin(beta2m));

 %- Variation of true longitude DL2(x)
 sin_DL2 = 1/sin_alpha2 * sin(delta2);
 cos_DL2 = cos(psi*(1-x)) .* cos(delta2);
 DL2 = atan2(sin_DL2,cos_DL2);
 
 sin_DL2m = 1/sin_alpha2 * sin(delta2m);
 cos_DL2m = cos(psi*(1-xm)) .* cos(delta2m);
 DL2m = atan2(sin_DL2m,cos_DL2m);
 
 %- Derivatives of DL2(x) wrt x : DL2_x DL2_xx DL2_xxx
 DL2_x    = delta2_x./(sin_alpha2*cos(psi*(1-x)));
 DL2_xx   = (delta2_xx - psi*sin_alpha2*sin(psi*(1-x)).*DL2_x)./(sin_alpha2*cos(psi*(1-x)));
 
 DL2m_x    = delta2m_x./(sin_alpha2*cos(psi*(1-xm)));
 DL2m_xx   = (delta2m_xx - psi*sin_alpha2*sin(psi*(1-xm)).*DL2m_x)./(sin_alpha2*cos(psi*(1-xm)));
          
          
%-------------------------------------------------------------------------%
% ATTRACTOR DISTANCE
 %- True long on the initial and final orbit: l1(x) l1_x(x) l2(x) l2_x(x)
%  l1 = L1 + DL1; 
%  l1m = L1 + DL1m; 
%  
%  l2 = L2 + DL2;
%  l2m = L2 + DL2m;
 
 cosl1 = cos(DL1)*cos(L1) - sin(DL1)*sin(L1);
 cosl2 = cos(DL2)*cos(L2) + sin(DL2)*sin(L2);
 
 cosl1m = cos(DL1m)*cos(L1) - sin(DL1m)*sin(L1);
 cosl2m = cos(DL2m)*cos(L2) + sin(DL2m)*sin(L2);
 
 sinl1 = sin(DL1)*cos(L1) + cos(DL1)*sin(L1);
 sinl2 = sin(DL1)*cos(L1) - cos(DL1)*sin(L1);
 
 sinl1m = sin(DL1m)*cos(L1) + cos(DL1m)*sin(L1);
 sinl2m = sin(DL1m)*cos(L1) - cos(DL1m)*sin(L1);
 
 
 %- Quantity qi(x) i = 1,2 needed to compute de distance from the attractor

 q1     = 1 + f1*cosl1 + g1*sinl1;
 q1_x   = (- f1*sinl1 + g1*cosl1).*DL1_x;
 q1_xx  = (1 - q1).*DL1_x.^2 + q1_x.*DL1_xx./DL1_x.^2;
 
 q1m     = 1 + f1*cosl1m + g1*sinl1m;
 q1m_x   = (- f1*sinl1m + g1*cosl1m).*DL1m_x;
 q1m_xx  = (1 - q1m).*DL1m_x.^2 + q1m_x.*DL1m_xx./DL1m_x.^2;

 
 q2     = 1 + f2*cosl2 + g2*sinl2;
 q2_x   = (- f2*sinl2 + g2*cosl2).*DL2_x;
 q2_xx  = (1 - q2).*DL2_x.^2 + q2_x.*DL2_xx./DL2_x.^2;

 q2m     = 1 + f2*cosl2m + g2*sinl2m;
 q2m_x   = (- f2*sinl2m + g2*cosl2m).*DL2m_x;
 q2m_xx  = (1 - q2m).*DL2m_x.^2 + q2m_x.*DL2m_xx./DL2m_x.^2;
 
 %- Distance s/c from the attractor (through MEE definition) : si(x) and
 %  derivatives si_x(x) si_xx(x) si_xxx(x) where i = 1,2
 
 s1     = p1./q1;
 s1_x   = -p1*q1_x./(q1.^2);
 s1_xx  = 2*p1*q1_x.^2./(q1.^3) - p1*q1_xx./(q1.^2);
 
 s1m     = p1./q1m;
 s1m_x   = -p1*q1m_x./(q1m.^2);
 s1m_xx  = 2*p1*q1m_x.^2./(q1m.^3) - p1*q1m_xx./(q1m.^2);
 
 s2     = p2./q2;
 s2_x   = -p2*q2_x./(q2.^2);
 s2_xx  = 2*p2*q2_x.^2./(q2.^3) - p2*q2_xx./(q2.^2);
 
 s2m     = p2./q2m;
 s2m_x   = -p2*q2m_x./(q2m.^2);
 s2m_xx  = 2*p2*q2m_x.^2./(q2m.^3) - p2*q2m_xx./(q2m.^2);
 
%-------------------------------------------------------------------------%
% Interpolating function 
xi     = a*x.^8 - (20 + 4*a)*x.^7 + (70 + 6*a)*x.^6 - (84 + 4*a)*x.^5 + (35 + a)*x.^4;
xi_x   = 8*a*x.^7 - 7*(20+4*a)*x.^6 + 6*(70 + 6*a)*x.^5 - 5*(84 + 4*a)*x.^4 + 4*(35 + a)*x.^3;
xi_xx  = 56*a*x.^6 - 42*(20+4*a)*x.^5 + 30*(70 + 6*a)*x.^4 - 20*(84 + 4*a)*x.^3 + 12*(35 + a)*x.^2;

xim     =  a*xm.^8 - (20 + 4*a)*xm.^7 + (70 + 6*a)*xm.^6 - (84 + 4*a)*xm.^5 + (35 + a)*xm.^4;
xim_x   =  8*a*xm.^7 - 7*(20+4*a)*xm.^6 + 6*(70 + 6*a)*xm.^5 - 5*(84 + 4*a)*xm.^4 + 4*(35 + a)*xm.^3;
xim_xx  =  56*a*xm.^6 - 42*(20+4*a)*xm.^5 + 30*(70 + 6*a)*xm.^4 - 20*(84 + 4*a)*xm.^3 + 12*(35 + a)*xm.^2;

% Ds Ds_x  Ds_xx 
Ds     = s2 - s1;
Ds_x   = s2_x - s1_x;
Ds_xx  = s2_xx - s1_xx;

Dsm     = s2m - s1m;
Dsm_x   = s2m_x - s1m_x;
Dsm_xx  = s2m_xx - s1m_xx;

% Ddelta Ddelta_x Ddelta_xx Ddelta_xxx 
Ddelta     = delta2 - delta1 ;
Ddelta_x   = delta2_x - delta1_x;
Ddelta_xx  = delta2_xx - delta1_xx;

Ddeltam     = delta2m - delta1m ;
Ddeltam_x   = delta2m_x - delta1m_x;
Ddeltam_xx  = delta2m_xx - delta1m_xx;

% Sun distance s(x) and declination delta(x) + derivatives
s     =  (s2 - s1) .* xi + s1;
delta =  (delta2 - delta1) .* xi + delta1; 

s_x   =  Ds_x.*xi + xi_x .*Ds + s1_x;
s_xx  =  Ds_xx.*xi + xi_xx.*Ds + 2*Ds_x.*xi_x + s1_xx;

delta_x   = Ddelta_x.*xi + xi_x.*Ddelta + delta1_x;
delta_xx  = Ddelta_xx.*xi + xi_xx.*Ddelta + 2*Ddelta_x.*xi_x + delta1_xx;

sm     = (s2m - s1m).* xim + s1m;
deltam = (delta2m - delta1m) .* xim + delta1m; 

sm_x   = Dsm_x.*xim + xim_x .*Dsm + s1m_x;
sm_xx  = Dsm_xx.*xim + xim_xx.*Dsm + 2*Dsm_x.*xim_x + s1m_xx;

deltam_x   = Ddeltam_x.*xim + xim_x.*Ddeltam + delta1m_x;
deltam_xx  = Ddeltam_xx.*xim + xim_xx.*Ddeltam + 2*Ddeltam_x.*xim_x + delta1m_xx;

% r and its derivatives 
r     = s.*cos(delta);
r_x   = - s.*delta_x.*sin(delta) + s_x .* cos(delta);
r_xx  = - (2*s_x.*delta_x + s_x.*delta_xx).*sin(delta) + (s_xx - s.*delta_x.^2).*cos(delta);

rm     = sm.*cos(deltam);
rm_x   = - sm.*deltam_x.*sin(deltam) + sm_x .* cos(deltam);
rm_xx  = - (2*sm_x.*deltam_x + sm_x.*deltam_xx).*sin(deltam) + (sm_xx - sm.*deltam_x.^2).*cos(deltam);

% Derivative of x
x_t_2 = sim.mu*r./(s.^3.*(r*psi^2 - r_xx + 2*r_x.^2./r));  
x_t   = sqrt(x_t_2);

xm_t_2 = abs(sim.mu*rm./(sm.^3.*(rm*psi^2 - rm_xx + 2*rm_x.^2./rm)));  
xm_t   = sqrt(xm_t_2);

% ----------------------------------------------------------------------- %
% dTOF at each theta 
dTOF  = 1./x_t;
dTOFm = 1./xm_t;

% step size
h = x(2)- x(1);

% Cavalieri-Simpson method
I = 0;
for i = 2:n_sol
    I = I + dTOF(i-1) + 4*dTOFm(i-1) + dTOF(i);
end

I = I*h/6;

% Computation of TOF error (residual)
error_TOF = abs(TOF-I)/abs(TOF) %% abs l'ho aggiunto io
err(j) = error_TOF;

end
end
