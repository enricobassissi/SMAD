function [output] = NL_interpolator( RI , RF , VI , VF , Nrev , TOF ,M ,Isp ,sim )

%[output ] = CW_LowLambert( RI , RF , VI , VF , Nrev , TOF ,M ,hp , kp , PS ,sim )

%-- Input:
   % RI RF VI VF
   % Nrev
   % TOF
   % M
   % Isp
   % sim.n_sol = 100;                    % number of computational nodes
   % sim.x = linspace(0,1,sim.n_sol)';   % non dimensional anomaly (varies in
                                         % [0,1])
   % sim.g0 
 
% ----------------------------------------------------------------------- % 

%-- Measurament units definition %%%%%%%%
 DU = sim.DU ;                    % Distance unit lenght[km]
 TU = sim.TU;                     % Time unit duration [s]

 x = sim.x;
 n_sol = sim.n_sol;
 
%-- Reference System 
 % magnitude of initial and final position vectors
 ri = sqrt(RI(1)*RI(1) + RI(2)*RI(2) + RI(3)*RI(3));
 rf = sqrt(RF(1)*RF(1) + RF(2)*RF(2) + RF(3)*RF(3));

 % initial and final position versors
 RIv = RI/ri;
 RFv = RF/rf;
 
 % vectorial product between RIv and RFv and its normal 
 RIv_cross_RFv      = [ RIv(2)*RFv(3) - RIv(3)*RFv(2);
                        RIv(3)*RFv(1) - RIv(1)*RFv(3);
                        RIv(1)*RFv(2) - RIv(2)*RFv(1)];
                       
 norm_RIv_cross_RFv = sqrt(RIv_cross_RFv(1)*RIv_cross_RFv(1) + RIv_cross_RFv(2)*RIv_cross_RFv(2) ...
                      + RIv_cross_RFv(3)*RIv_cross_RFv(3));
                       
 % scalar product between RIv and RFv
 RIv_dot_RFv = RIv(1)*RFv(1) + RIv(2)*RFv(2) + RIv(3)*RFv(3);
 
 % initial and final angular momentum vector direction
 HI = cross(RI,VI)/norm(cross(RI,VI));
 HF = cross(RF,VF)/norm(cross(RF,VF));
 
 % Transfer plane definition (z = 0) [direction perpendicular to RI, RF]
 % and transfer angle psi
 Href = RIv_cross_RFv/ norm_RIv_cross_RFv;
 
 
    % Transfer must be counterclockwise: 
    if Href(3) < 0 %%% paper usa condizione diversa : dot(RIv_cross_RFv, HI) < 0 
        Href = - Href;
    end
    
    % Singular case:  
    if  norm_RIv_cross_RFv < 1e-10 %%% paper usa condizione diversa : dot(RIv_cross_RFv, HI) = 0 
        Href = (HI + HF)/norm(HI + HF); %%%% prof la scrive diversa (considera solo HI)
        
        if RIv_dot_RFv > 0
            psi1 = 0;
        else
            psi1 = pi;
        end
    
    else
        if RIv_cross_RFv(3) > 0   %%% paper usa condizione diversa :  cross(RIv_cross_RFv, Href) > 0
            psi1 = acos(RIv_dot_RFv);
        elseif RIv_cross_RFv(3) < 0  %%% paper usa condizione diversa :  cross(RIv_cross_RFv, Href) < 0
            psi1 = 2*pi - acos(RIv_dot_RFv);
        end
    end

    
 % For more than 1 revolution:
 psi = psi1 + 2*pi*Nrev;   
 

%-- theta(t) vector 
theta = psi*sim.x;  

% In ConWay: x is defined outside (vector dimension n_sol from 0 to 1)
% Inside he defines l = L*sim.x where L = psi

% in the paper: theta given and inside x = theta/psi;


% ----------------------------------------------------------------------- %
%-- Initial modified equinoctial elements MEE1 = [p1 f1 g1 h1 k1 L1]
 [coe1] = coe_from_sv(RI,VI,sim.mu);
 a1 = coe1(1); e1 = coe1(2); i1 = coe1(3); OM1 = coe1(4); om1 = coe1(5); th1 = coe1(6);

 p1 = a1*(1-e1^2);
 f1 = e1*cos(om1) + OM1;
 g1 = e1*sin(om1) + OM1;
 h1 = tanh(i1/2) * cos(OM1);
 k1 = tanh(i1/2) * sin(OM1);
 L1 = th1 + om1 + OM1;

%-- Final modified equinoctial elements MEE2 = [p2 f2 g2 h2 k2 L2]
 [coe2] = coe_from_sv(RF,VF,sim.mu)
 a2 = coe2(1); e2 = coe2(2); i2 = coe2(3); OM2 = coe2(4); om2 = coe2(5); th2 = coe2(6);

 p2 = a2*(1-e2^2);
 f2 = e2*cos(om2) + OM2;
 g2 = e2*sin(om2) + OM2;
 h2 = tanh(i2/2) * cos(OM2);
 k2 = tanh(i2/2) * sin(OM2);
 L2 = th2 + om2 + OM2;

% ----------------------------------------------------------------------- %
%-- Departure orbit

 %- Inclination of the departure orbit above the reference frame (alpha1)
 cos_alpha1 = dot(HI,Href); % dot product 
 
 if dot(VI,Href) > 0 %%% dot product - caso uguale a zero? vi vf sono velocità iniziali e finali? %%%%%%%
     csi1 = 1;
 elseif dot(VI,Href) < 0
     csi1 = -1;
 end
 sin_alpha1 = csi1*sqrt(1 - cos_alpha1^2);
 
 alpha1 = atan(sin_alpha1/cos_alpha1);
 
 % CONTROLLA I .* E ./
 
 %- Angle beta1(x) and derivatives needed to compute the declination
 beta1     = acos(sin_alpha1*cos(psi*x));
 beta1_x   = psi*sin_alpha1*sin(psi*x)./sin(beta1);
 beta1_xx  = (psi^2*sin_alpha1*cos(psi*x) - cos(beta1).*beta1_x.^2)./sin(beta1);
 beta1_xxx = (-psi^3*sin_alpha1*sin(psi*x) - 3*beta1_x.*beta1_xx.*cos(beta1) + beta1_x.^3.*sin(beta1))./sin(beta1);
 
 %- Declination above the ref. plane delta1(x) - è univoc definita?? %%%%
 delta1 = asin(sin_alpha1*sin(psi*x)./sin(beta1));
 
 %- Derivatives of delta1(x) wrt x : delta1_x delta1_xx delta1_xxx
 delta1_x   = (psi *sin_alpha1 * cos(psi*x) - beta1_x.*cos(beta1).*sin(delta1))./(cos(delta1) .* sin(beta1));
 delta1_xx  = (- psi^2 *sin_alpha1 .* sin(psi*x) ...
              +  sin(beta1).*sin(delta1).*(delta1_x.^2 + beta1_x.^2))./(cos(delta1) .* sin(beta1)) ...
              - (2*delta1_x.*beta1_x.*cos(delta1).*cos(beta1) ...
              +  beta1_xx.*sin(delta1).*cos(beta1))./(cos(delta1) .* sin(beta1));
 delta1_xxx = (- psi^3*sin_alpha1*cos(psi*x) ... 
              +  sin(beta1).*sin(delta1).*(3*delta1_x.*delta1_xx + 3*beta1_x.*beta1_xx))./(cos(delta1).*sin(beta1)) ...
              +  (- cos(beta1).*cos(delta1).*(3*delta1_x.*delta1_xx + 3*beta1_x.*beta1_xx) ...
              + cos(beta1).*sin(delta1).*(3*delta1_x.^2 .*beta1_x + beta1_x.^3 - beta1_xxx))./(cos(delta1).*sin(beta1)) ...
              + sin(beta1).*cos(delta1).*(3*beta1_x.^2 .*delta1_x + delta1_x.^3)./(cos(delta1).*sin(beta1));
 
 %- Variation of true longitude DL1(x)
 sin_DL1 = 1/sin_alpha1 * sin(delta1);
 cos_DL1 = cos(psi*x) .* cos(delta1);
 DL1 = atan(sin_DL1./cos_DL1);
 
 %- Derivatives of DL1(x) wrt x * DL1_x DL1_xx DL1_xxx
 DL1_x    = delta1_x./(sin_alpha1.*cos(psi*x));
 DL1_xx   = (delta1_xx + psi*sin_alpha1.*sin(psi*x) .*DL1_x)./(sin_alpha1.*cos(psi*x));
 DL1_xxx  = (delta1_xxx + 2*psi*sin_alpha1*sin(psi*x).*DL1_xx + psi^2*sin_alpha1*cos(psi*x).*DL1_x)./(sin_alpha1.*cos(psi*x));
 

 
%-- Arrival orbit

 %- Inclination of the arrival orbit above the reference frame (alpha1)
 cos_alpha2 = dot(HF,Href); % dot product 
 
 if dot(VF,Href) < 0 %%% dot product - uguale a zero?
     csi2 = 1;
 elseif dot(VF,Href) > 0
     csi2 = -1;
 end
 sin_alpha2 = csi2*sqrt(1 - cos_alpha2^2);
 
 alpha2 = atan(sin_alpha2/cos_alpha2);
 
 %- Angle beta2(x) and derivatives needed to compute the declination
 beta2     = acos(sin_alpha2*cos(psi*(1 -x)));
 beta2_x   = - psi*sin_alpha2*sin(psi*(1-x))./sin(beta2);
 beta2_xx  = (- psi^2*sin_alpha2.*cos(psi*(1-x)) - cos(beta2).*beta2_x.^2)./sin(beta2);
 beta2_xxx = (psi^3*sin_alpha2*sin(psi*(1-x)) - 3*beta2_x.*beta2_xx.*cos(beta2) + beta2_x.^3.*sin(beta2))./sin(beta2);
 
 %- Declination above the ref. plane delta1(x) - è univoc definita?? %%%%
 delta2 = asin(sin_alpha2*sin(psi*(1-x))./sin(beta2));
 
 %- Derivatives of delta2(x) wrt x : delta2_x delta2_xx delta2_xxx
 delta2_x   = (- psi *sin_alpha2 * cos(psi*(1-x)) - beta2_x.*cos(beta2).*sin(delta2))./(cos(delta2) .* sin(beta2));
 delta2_xx  = (- psi^2 *sin_alpha2 * sin(psi*(1-x)) ...
              +  sin(beta2).*sin(delta2).*(delta2_x.^2 + beta2_x.^2))./(cos(delta2) .* sin(beta2)) ...
              - (2*delta2_x.*beta2_x .*cos(delta2) .*cos(beta2) ...
              +  beta2_xx .*sin(delta2) .*cos(beta2))./(cos(delta2) .* sin(beta2));
 delta2_xxx = (  psi^3*sin_alpha2*cos(psi*(1-x)) ... 
              +  sin(beta2).*sin(delta2).*(3*delta2_x.*delta2_xx + 3*beta2_x.*beta2_xx))./(cos(delta2).*sin(beta2)) ...
              +  (- cos(beta2).*cos(delta2).*(3*delta2_x .*delta2_xx + 3*beta2_x.*beta2_xx) ...
              + cos(beta2).*sin(delta2).*(3*delta2_x.^2.*beta2_x + beta2_x.^3 - beta2_xxx))./(cos(delta2).*sin(beta2)) ...
              + sin(beta2).*cos(delta2).*(3*beta2_x.^2.*delta2_x + delta2_x.^3)./(cos(delta2).*sin(beta2));
 
 %- Variation of true longitude DL2(x)
 sin_DL2 = 1/sin_alpha2 * sin(delta2);
 cos_DL2 = cos(psi*(1-x)) .* cos(delta2);
 DL2 = atan(sin_DL2./cos_DL2);
 
 %- Derivatives of DL2(x) wrt x : DL2_x DL2_xx DL2_xxx
 DL2_x    = delta2_x./(sin_alpha2*cos(psi*(1-x)));
 DL2_xx   = (delta2_xx - psi*sin_alpha2*sin(psi*(1-x)).*DL2_x)./(sin_alpha2*cos(psi*(1-x)));
 DL2_xxx  = (delta2_xxx + 2*psi*sin_alpha2*sin(psi*(1-x)).*DL2_xx ...
            + psi^2*sin_alpha2*cos(psi*(1-x)).*DL2_x)./(sin_alpha2*cos(psi*(1-x))); 
 

 
%-- Attractor distance
 %- True long on the initial and final orbit: l1(x) l1_x(x) l2(x) l2_x(x)
 l1 = L1 + DL1; %%%% L1 è 1x1 e DL1 è 100x1 ?????
 l1_x = DL1_x;
 
 l2 = L2 + DL2; %%%% L2 è 1x1 e DL2 è 100x1 ?????
 l2_x = -DL2_x;
 
 %- Quantity qi(x) i = 1,2 needed to compute de distance from the attractor

 q1     = 1 + f1*cos(l1) + g1*sin(l1);
 q1_x   = (- f1*sin(l1) + g1*cos(l1)).*DL1_x;
 q1_xx  = (1 - q1).*DL1_x.^2 + q1_x.*DL1_xx./DL1_x.^2;
 q1_xxx = - q1_x.*DL1_x.^2 + 3*(1-q1).*DL1_x.*DL1_xx + q1_x.*DL1_xxx./DL1_x;
 
 q2     = 1 + f2*cos(l2) + g1*sin(l2);
 q2_x   = (- f2*sin(l2) + g2*cos(l2)).*DL2_x;
 q2_xx  = (1 - q2).*DL2_x.^2 + q2_x.*DL2_xx./DL2_x.^2;
 q2_xxx = - q2_x.*DL2_x.^2 + 3*(1-q2).*DL2_x.*DL2_xx + q2_x.*DL2_xxx./DL2_x;
 
 %- Distance s/c from the attractor (through MEE definition) : si(x) and
 %  derivatives si_x(x) si_xx(x) si_xxx(x) where i = 1,2
 
 s1     = p1./q1;
 s1_x   = -p1*q1_x./(q1.^2);
 s1_xx  = 2*p1*q1_x.^2./(q1.^3) - p1*q1_xx./(q1.^2);
 s1_xxx = -6*p1*q1_x.^3./(q1.^4) + 6*p1*q1_x.*q1_xx./(q1.^3) - p1.*q1_xxx./(q1.^2);
 
 s2     = p2./q2;
 s2_x   = -p2*q2_x./(q2.^2);
 s2_xx  = 2*p2*q2_x.^2./(q2.^3) - p2*q2_xx./(q2.^2);
 s2_xxx = -6*p2*q2_x.^3./(q2.^4) + 6*p2*q2_x.*q2_xx./(q2.^3) - p2*q2_xxx./(q2.^2);
 
% ----------------------------------------------------------------------- %
%- Interpolationg function coefficient a - fsolve + Cavalieri-Simpson to
%  solve the integral
a0  = 0;
% fun = @(a) find_a(x,psi,TOF, s1,s1_x, s1_xx, s2,s2_x,s2_xx, delta1,delta1_x,delta1_xx,delta2,delta2_x,delta2_xx)
% options=optimoptions('fsolve', 'TolFun', 1e-13, 'TolX', 1e-13,'Display','off');
% a = fsolve(fun,a0,options); %% alternative: fzero 
a = 0;

%- Interpolating function xi(x,a) and derivatives
 % Initial and final position and velocity shall be constrained, therefore
 % the function shall be continuous (between 0 and 1).
 % If you want TOF free (not imposed) : a = 0.      %%%% controlla derivate
xi     = a*x.^8 - (20 + 4*a)*x.^7 + (70 + 6*a)*x.^6 - (84 + 4*a)*x.^5 + (35 + a)*x.^4;
xi_x   = 8*a*x.^7 - 7*(20+4*a)*x.^6 + 6*(70 + 6*a)*x.^5 - 5*(84 + 4*a)*x.^4 + 4*(35 + a)*x.^3;
xi_xx  = 56*a*x.^6 - 42*(20+4*a)*x.^5 + 30*(70 + 6*a)*x.^4 - 20*(84 + 4*a)*x.^3 + 12*(35 + a)*x.^2;
xi_xxx =  336*a*x.^5 - 210*(20+4*a)*x.^4 + 120*(70 + 6*a)*x.^3 - 60*(84 + 4*a)*x.^2 + 24*(35 + a)*x;

%- Sun distance s(x) and declination delta(x) 
s     = (s2 - s1).* xi + s1;
delta = (delta2 - delta1).* xi + delta1; 


%-- Transformation between spherical and cylindrical coordinates
r = s.*cos(delta);
z = s.*sin(delta);

% Ds Ds_x  Ds_xx Ds_xxx %%%% controlla
Ds     = s2 - s1;
Ds_x   = s2_x - s1_x;
Ds_xx  = s2_xx - s1_xx;
Ds_xxx = s2_xxx - s1_xxx;


% Ddelta Ddelta_x Ddelta_xx Ddelta_xxx %%%% controlla
Ddelta     = delta2 - delta1 ;
Ddelta_x   = delta2_x - delta1_x;
Ddelta_xx  = delta2_xx - delta1_xx;
Ddelta_xxx = delta2_xxx - delta1_xxx;


% ----------------------------------------------------------------------- %
%-- Derivatives wrt x needed  for the parametrized dynamical system

% s_x s_xx s_xxx
s_x   = Ds_x.*xi + xi_x .*Ds + s1_x;
s_xx  = Ds_xx.*xi + xi_xx.*Ds + 2*Ds_x.*xi_x + s1_xx;
s_xxx = Ds_xxx.*xi + xi_xxx.*Ds + 3*Ds_xx.*xi_x + 3*Ds_x.*xi_xx + s1_xxx;

% delta_x delta_xx delta_xxx
delta_x   = Ddelta_x.*xi + xi_x.*Ddelta + delta1_x;
delta_xx  = Ddelta_xx.*xi + xi_xx.*Ddelta + 2*Ddelta_x.*xi_x + delta1_xx;
delta_xxx = Ddelta_xxx.*xi + xi_xxx.*Ddelta +3*Ddelta_xx.*xi_x + 3*Ddelta_x.*xi_xx + delta1_xxx;

% r_x r_xx r_xxx
r_x   = - s.*delta_x + s_x .* cos(delta);
r_xx  = - (2*s_x.*delta_x + s_x.*delta_xx).*sin(delta) + (s_xx - s.*delta_x.^2).*cos(delta);
r_xxx = - (3*(s_xx.*delta_x + s_x.*delta_xx) + s.*(delta_xxx - delta_x.^3)).*sin(delta) + ...
          (s_xxx - 3*delta_x.*(s_x.*delta_x + s.*delta_xx)).*cos(delta);

% z_x z_xx
z_x   = s.*delta_x.*cos(delta) + s_x.*sin(delta);
z_xx  = (2*s_x.*delta_x + s.*delta_xx).*cos(delta) + (s_xx -s.*delta_x.^2).*sin(delta);

% ----------------------------------------------------------------------- %

%-- Parametrized EOM

%- Square of derivative of x wrt time : x_t_2
x_t_2 = abs(sim.mu*r./(s.^3.*(r*psi^2 - r_xx + 2*r_x.^2./r)));  %%% CONTROLLA -> abs l'ho aggiunto io
x_t   = sqrt(x_t_2);

%- Quantities to compute the second derivative of x wrt time
Nu = sim.mu*r; %%% controlla
De = s.^3.*(r*psi^2 - r_xx + 2*r_x.^2./r); %%% controlla

Nu_x = sim.mu*r_x;
De_x = 3*s_x./s.*De + s.^3.*(r_x*psi^2 - r_xxx + (2*r.*r_x.*r_xx - r_x.^3)./(r.^2));

%- Second derivative of x wrt time: x_tt
x_tt = 0.5*(Nu_x - x_t_2.*De_x)./De;
   
% ----------------------------------------------------------------------- %

%- Thrust angle gamma -> It is imposed equal to the flight path angle in
%  order to have the in-plane thrust imposed as tangential only. In this way
%  the x time derivative can be analytically computed removing the dependency
%  from thrust per unit mass. 
gamma = atan(r_x./(r*psi));

%- In plane thrust per unit mass
Tin2m  = 1/cos(gamma) * (2*psi*r_x.*x_t_2 + r.*psi.*x_tt);

%- Out of plane thrust per unit mass
Tout2m = z_xx.*x_t_2 + z_x.*x_tt + sim.mu./(s.^3).*z;

%- Thrust (absolute value) per unit mass
T2m = sqrt(Tin2m.^2 + Tout2m.^2);

%- Mass time history from Tsiolkovsky equation: m_t
   % I consider the backward case (from dry to wet) %%%%%%%%%%

dtheta = theta(2)- theta(1);
dtheta2 = dtheta/2;
dtheta12 = dtheta/12;
dtheta24 = dtheta/24;

   if sim.direction == -1   
       m = theta; %% perchè?
       m(sim.n_sol) = M;
       K = - T2m./(Isp*sim.g0);
       m(sim.n_sol-1) = m(sim.n_sol) -   dtheta2*(K(sim.n_sol)*m(sim.n_sol)...
                        +K(sim.n_sol-1)*m(sim.n_sol-1)); % (?) Alternative: FE
       m(sim.n_sol-2) = m(sim.n_sol-1) - dtheta2*(K(sim.n_sol-1)*m(sim.n_sol-1)...
                        +K(sim.n_sol-2)*m(sim.n_sol-2));
        
       for i = sim.n_sol-2:-1:2 % predictor corrector AB3AM4 
           m(i-1) = m(i) - dtheta12*(23*K(i)*m(i) - 16*K(i+1)*m(i+1) + 5*K(i+2)*m(i+2) ) ;
           m(i-1) = m(i) - dtheta24*(9*K(i-1)*m(i-1)+19*K(i)*m(i) - 5*K(i+1)*m(i+1) + K(i+2)*m(i+2) ) ;
       end
   end

   
%- time vector
d_time =1./x_t;
dtheta6 = dtheta./6;
t(1) = 0; % Initialization

for i = 2:n_sol-1
    t(i) = t(i-1) + dtheta6 * (d_time(i-1) + 4*d_time(i) + d_time(i+1) ) ; %% RK4 (?)
end
t(n_sol) = t(n_sol-1) + d_time(n_sol-1)*dtheta ;
    
output.m       = m ;
output.t       = t ;    

end