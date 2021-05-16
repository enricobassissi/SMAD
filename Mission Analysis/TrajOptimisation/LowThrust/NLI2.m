function [output] = NLI2( RI , RF , VI , VF , N_REV , TOF ,M ,Isp ,sim)
%% Input:
   % RI RF VI VF
   % Nrev
   % TOF
   % M
   % Isp
   % sim.n_sol = 100;                    % number of computational nodes
   % sim.x = linspace(0,1,sim.n_sol)';   % non dimensional anomaly (varies in
                                         % [0,1])
   % sim.g0 
   % sim.DU sim.TU
   % sim.mu_dim                          % dimensional mu
 
%% Measurament units definition 
 DU = sim.DU ;                    % Distance unit lenght[km]
 TU = sim.TU;                     % Time unit duration [s]

 x = sim.x;
 n_sol = sim.n_sol;
 
%% MEE of the orbits
%-- Initial modified equinoctial elements MEE1 = [p1 f1 g1 h1 k1 L1]
 coe1 = coe_from_sv(RI*DU,VI*DU/TU,sim.mu_dim);  
 a1 = coe1(1); e1 = coe1(2); i1 = coe1(3); OM1 = coe1(4); om1 = coe1(5); th1 = coe1(6);

 p1 = a1*(1-e1^2)/DU;
 f1 = e1*cos(om1 + OM1); 
 g1 = e1*sin(om1 + OM1); 
 h1 = tan(i1/2) * cos(OM1);
 k1 = tan(i1/2) * sin(OM1);
 L1 = th1 + om1 + OM1; 
 
 sinL1 = sin(L1);
 cosL1 = cos(L1);
 
%-- Final modified equinoctial elements MEE2 = [p2 f2 g2 h2 k2 L2]
 coe2 = coe_from_sv(RF*DU,VF*DU/TU, sim.mu_dim);
 a2 = coe2(1); e2 = coe2(2); i2 = coe2(3); OM2 = coe2(4); om2 = coe2(5); th2 = coe2(6);

 p2 = a2*(1-e2^2)/DU;
 f2 = e2*cos(om2 + OM2); 
 g2 = e2*sin(om2 + OM2); 
 h2 = tan(i2/2) * cos(OM2);
 k2 = tan(i2/2) * sin(OM2);
 L2 = th2 + om2 + OM2;
 
 sinL2 = sin(L2);
 cosL2 = cos(L2);
 
%% General part
 % initial and final 3d positions and velocities
 qi = 1 + f1*cosL1 + g1*sinL1;
 ri = p1/qi;
 alpha2i = h1^2 - k1^2;
 chi2i = h1^2 + k1^2;
 s2i = 1 + chi2i;
 due_hki = 2*h1*k1;
 
 RI = ri/s2i*[ cosL1*(1+alpha2i) + due_hki*sinL1;
               sinL1*(1-alpha2i) + due_hki*cosL1;
               2*(h1*sinL1 - k1*cosL1) ];
   
 VI = (sim.mu/p1)^0.5/s2i * [ -(sinL1*(1+alpha2i) - due_hki*cosL1 + g1 - f1*due_hki + alpha2i*g1);
            -(cosL1*(alpha2i-1)+due_hki*sinL1 - f1 + g1*due_hki + alpha2i*f1);
            2*(h1*cosL1 + k1*sinL1 + f1*h1 + g1*k1)];
 
 vi = sqrt(VI(1)^2+VI(2)^2+VI(3)^2);

 qf = 1 + f2*cosL2 + g2*sinL2;
 rf = p2/qf;
 alpha2f = h2^2 - k2^2;
 chi2f = h2^2 + k2^2;
 s2f = 1 + chi2f;
 due_hkf = 2*h2*k2;
 
 RF = rf/s2f*[ cosL2*(1+alpha2f) + due_hkf*sinL2;
               sinL2*(1-alpha2f) + due_hkf*cosL2;
               2*(h2*sinL2 - k2*cosL2) ];
   
 VF = (sim.mu/p2)^0.5/s2f * [ -(sinL2*(1+alpha2f) - due_hkf*cosL2 + g2 - f2*due_hkf + alpha2f*g2);
            -(cosL2*(alpha2f-1)+due_hkf*sinL2 - f2 + g2*due_hkf + alpha2f*f2);
            2*(h2*cosL2 + k2*sinL2 + f2*h2 + g2*k2)];
            
 vf = sqrt(VF(1)^2+VF(2)^2+VF(3)^2);
 
 % initial and final position versors
 RI_vers = RI/ri;
 RF_vers = RF/rf;
 
 VI_vers = VI/vi;
 VF_vers = VF/vf;
 
 % initial and final angular momentum and directions
 HI = [ RI(2)*VI(3) - RI(3)*VI(2);
        RI(3)*VI(1) - RI(1)*VI(3);
        RI(1)*VI(2) - RI(2)*VI(1)];
    
 HF = [ RF(2)*VF(3) - RF(3)*VF(2);
        RF(3)*VF(1) - RF(1)*VF(3);
        RF(1)*VF(2) - RF(2)*VF(1)];
 
 HI_vers = HI/(sqrt(HI(1)^2 + HI(2)^2 + HI(3)^2));
 HF_vers = HF/(sqrt(HF(1)^2 + HF(2)^2 + HF(3)^2));
 
 % normal of reference plane
 RI_vers_X_RF_vers = [ RI_vers(2)*RF_vers(3) - RI_vers(3)*RF_vers(2);
                       RI_vers(3)*RF_vers(1) - RI_vers(1)*RF_vers(3);
                       RI_vers(1)*RF_vers(2) - RI_vers(2)*RF_vers(1)];
                       
 nrm_ref_plane = sqrt(RI_vers_X_RF_vers(1)^2 + RI_vers_X_RF_vers(2)^2 + RI_vers_X_RF_vers(3)^2);
 if nrm_ref_plane > sim.tol_vers
     h_ref_vers = RI_vers_X_RF_vers/nrm_ref_plane;
 else
     h_ref_vers = 0.5*(HI_vers + HF_vers);
     h_ref_vers = h_ref_vers/(sqrt(h_ref_vers(1)^2+h_ref_vers(2)^2+h_ref_vers(3)^2));
 end
 
 if h_ref_vers(1)*HI_vers(1)+h_ref_vers(2)*HI_vers(2)+h_ref_vers(3)*HI_vers(3) < 0 % ccw rotation
     h_ref_vers = -h_ref_vers;
 end
 
 % angle (scalar product) between RIv and RFv
 RI_vers_DOT_RF_vers = RI_vers(1)*RF_vers(1) + RI_vers(2)*RF_vers(2) + RI_vers(3)*RF_vers(3);
 
 if sqrt(RI_vers_X_RF_vers(1)^2+RI_vers_X_RF_vers(2)^2+RI_vers_X_RF_vers(3)^2) < sim.tol_vers
     if RI_vers_DOT_RF_vers > 0
         PSY = 2*pi;
     else
         PSY = pi;
     end
 else
     if RI_vers_X_RF_vers(3) > 0
         PSY = acos(RI_vers_DOT_RF_vers) ; %+ 2*pi
     else
         PSY = 2*pi - acos(RI_vers_DOT_RF_vers);
     end
 end

 PSY = PSY + N_REV*2*pi;
 Px = PSY*x;
 PSY_2 = PSY^2;
 PSY_3 = PSY^3;
 
 % sine and cosine function of x
 sinP = sin(PSY);
 cosP = cos(PSY);
 
 sinPx = sin(Px);
 cosPx = cos(Px);
 
 sinPumx = sinP*cosPx - cosP*sinPx;
 cosPumx = sinP*sinPx + cosP*cosPx;
 
%% Parametrization of declination
 % a1 = Inclination of the departure orbit above the reference frame
 % initial orbit 
 cos_a1 = (h_ref_vers(1)*HI_vers(1) + h_ref_vers(2)*HI_vers(2) +  h_ref_vers(3)*HI_vers(3));
 sin_a1 = sqrt(1-cos_a1*cos_a1);
 
 if cos_a1 >=1 || cos_a1 <= -1
     d1 = 0*x;
     d1_x = d1;
     d1_xx = d1;
     d1_xxx = d1;
     cos_b1 = d1;
     b1_x = d1;
     b1_xx = d1;
     b1_xxx = d1;
     cosDL1 = cosPx;
     sinDL1 = sinPx;
     DL1_x = PSY;
     DL1_xx = 0;
     DL1_xxx = 0;
 else
     if (h_ref_vers(1)*VI_vers(1) + h_ref_vers(2)*VI_vers(2) + h_ref_vers(3)*VI_vers(3)) >= 0
         up1 = sin_a1; % departure orbit above ref plane
     else
         up1 = -sin_a1; % dep orbit below ref plane
     end
     cos_b1 = sin_a1*cosPx;
     sin_b1 = sqrt(1-cos_b1.*cos_b1);
     
     sin_d1 = up1*sinPx./sin_b1;
     d1 = asin(sin_d1);
     cos_d1 = cos(d1);
     
     sD1sB1 = sin_d1.*sin_b1;
     cD1sB1 = cos_d1.*sin_b1;
     sD1cB1 = sin_d1.*cos_b1;
     cD1cB1 = cos_d1.*cos_b1;
     
     b1_x = sin_a1*PSY*sinPx./sin_b1;
     b1_x_2 = b1_x.^2;
     b1_x_3 = b1_x.^3;
     b1_xx = (sin_a1*PSY_2*cosPx - b1_x_2.*cos_b1)./sin_b1;
     b1_xxx = (-PSY_3.*sin_a1.*sinPx - 3*b1_x.*b1_xx.*cos_b1 + b1_x_3.*sin_b1)./sin_b1;
     
     d1_x   = (up1*PSY*cosPx - b1_x.*sD1cB1)./(cD1sB1);
     d1_x_2 = d1_x.^2;
     d1_x_3 = d1_x.^3;
     d1_xx  = (- up1*PSY_2*sinPx + sD1sB1.*(b1_x_2 + d1_x_2) + ...
               - 2*d1_x.*b1_x.*cD1cB1 - b1_xx.*sD1cB1)./cD1sB1;
     d1_xxx = (- up1*PSY_3*cosPx + (3*sD1sB1.*(d1_x.*d1_xx + b1_x.*b1_xx)) +...
               - 3*cD1cB1.*(d1_x.*d1_xx + b1_x.*b1_xx) + ...
               sD1cB1.*(3*d1_x_2 .*b1_x + b1_x_3 - b1_xxx) + ...
               cD1sB1.*(3*b1_x_2 .*d1_x + d1_x_3))./cD1sB1;
    
     sinDL1 = sin_d1/up1;
     cosDL1 = cosPx.*cos_d1;

     denom1 = up1.*cosPx;
     DL1_x = d1_x./denom1; 
     DL1_xx = (d1_xx + PSY*up1.*sinPx.*DL1_x)./denom1;
     DL1_xxx = (d1_xxx + 2*PSY*up1*sinPx.*DL1_xx + PSY*up1*cosPx.*DL1_x)./denom1;
     
 end
 
% final orbit
cos_a2 = (h_ref_vers(1)*HF_vers(1) + h_ref_vers(2)*HF_vers(2) +  h_ref_vers(3)*HF_vers(3));
sin_a2 = sqrt(1-cos_a2*cos_a2);

if cos_a2 >=1 || cos_a2 <= -1
    d2 = 0*x;
    d2_x = d2;
    d2_xx = d2;
    d2_xxx = d2;
    cos_b2 = d2;
    b2_x = d2;
    b2_xx = d2;
    b2_xxx = d2;
    cosDL2 = cosPumx;
    sinDL2 = sinPumx;
    DL2_x = -PSY;
    DL2_xx = 0;
    DL2_xxx = 0;
else
    if (h_ref_vers(1)*VF_vers(1) + h_ref_vers(2)*VF_vers(2) + h_ref_vers(3)*VF_vers(3)) >= 0
        up2 = -sin_a2; % arr orbit below ref plane
    else
        up2 = sin_a2; % arr orbit above ref plane
    end
    cos_b2 = sin_a2*cosPumx;
    sin_b2 = sqrt(1-cos_b2.*cos_b2);

    sin_d2 = up2*sinPumx./sin_b2;
    d2 = asin(sin_d2);
    cos_d2 = cos(d2);

    sD2sB2 = sin_d2.*sin_b2;
    cD2sB2 = cos_d2.*sin_b2;
    sD2cB2 = sin_d2.*cos_b2;
    cD2cB2 = cos_d2.*cos_b2;

    b2_x = -(sin_a2*PSY*sinPumx./sin_b2);
    b2_x_2 = b2_x.^2;
    b2_x_3 = b2_x.^3;
    b2_xx = (sin_a2*PSY_2*cosPumx - b2_x_2.*cos_b2)./sin_b2;
    b2_xxx = (PSY_3.*sin_a2.*sinPumx - 3*b2_x.*b2_xx.*cos_b2 + b2_x_3.*sin_b2)./sin_b2;

    d2_x   = (-up2*PSY*cosPumx - b2_x.*sD2cB2)./(cD2sB2);
    d2_x_2 = d2_x.^2;
    d2_xx  = (- up2*PSY_2*sinPumx + sD2sB2.*(b2_x_2 + d2_x_2) + ...
           - 2*d2_x.*b2_x.*cD2cB2 - b2_xx.*sD2cB2)./cD2sB2;
    d2_xxx = (up2*PSY_3*cosPumx + (3*sD2sB2.*(d2_x.*d2_xx + b2_x.*b2_xx)) +...
           - 3*cD2cB2.*(d2_x.*d2_xx + b2_x.*b2_xx) + ...
           sD2cB2.*(3*d2_x_2 .*b2_x + b2_x_3 - b2_xxx) + ...
           cD2sB2.*(3*b2_x_2 .*d2_x + d1_x_3))./cD2sB2;

    sinDL2 = sin_d2/up2;
    cosDL2 = cosPumx.*cos_d2;
    
    denom2 = up2.*cosPumx;
    DL2_x = d2_x./denom2;
    DL2_xx = (d2_xx - PSY*up2.*sinPumx.*DL2_x)./denom2;
    DL2_xxx = (d2_xxx - 2*PSY*up2*sinPumx.*DL2_xx + PSY*up2*cosPumx.*DL2_x)./denom2;

end

Dd = d2 - d1;
Dd_x = d2_x - d1_x;
Dd_xx = d2_xx - d1_xx;
Dd_xxx = d2_xxx - d1_xxx;

cosL11 = cosDL1*cosL1 - sinDL1*sinL1; % L11 = DL1 + L1;
sinL11 = sinDL1*cosL1 + cosDL1*sinL1;

cosL22 = cosL2*cosDL2 + sinL2*sinDL2; % L22 = L2 - DL2;
sinL22 = sinL2*cosDL2 + cosL2*sinDL2;

%% Parametrisation of Attractor distance
% True long on the initial and final orbit: l1(x) l1_x(x) l2(x) l2_x(x)
%  l1 = L1 + DL1; 
%%%% L1 è 1x1 e DL1 è 100x1
%%%% L2 è 1x1 e DL2 è 100x1
dL1 = DL1_x;
dL1_2 = dL1.^2;
ddL1 = DL1_xx;
dddL1 = DL1_xxx;

dL2 = -DL2_x;
dL2_2 = dL2.^2;
ddL2 = -DL2_xx;
dddL2 = -DL2_xxx;

%- Quantity qi(x) i = 1,2 needed to compute de distance from the attractor
q1     = 1 + f1*cosL11 + g1*sinL11;
q1_x   = dL1.*(- f1*sinL11 + g1*cosL11);
q1_xx  = dL1_2.*(1 - q1) + q1_x.*ddL1./dL1;
q1_xxx = - q1_x.*dL1_2 + 3*(1-q1).*dL1.*ddL1 + q1_x.*dddL1./dL1;

q2     = 1 + f2*cosL22 + g2*sinL22;
q2_x   = dL2.*(- f2*sinL22 + g2*cosL22);
q2_xx  = dL2_2.*(1 - q2) + q2_x.*ddL2./dL2;
q2_xxx = - q2_x.*dL2_2 + 3*(1-q2).*dL2.*ddL2 + q2_x.*dddL2./dL2;

%- Distance s/c from the attractor (through MEE definition) : si(x) and
%  derivatives si_x(x) si_xx(x) si_xxx(x) where i = 1,2
s1     = p1./q1;
s1_x   = -p1*q1_x./(q1.^2);
s1_xx  = p1*(2*q1_x.^2./(q1.^3) - q1_xx./(q1.^2));
s1_xxx = p1*(-6*q1_x.^3./(q1.^4) + 6*q1_x.*q1_xx./(q1.^3) - q1_xxx./(q1.^2));

s2     = p2./q2;
s2_x   = -p2*q2_x./(q2.^2);
s2_xx  = p2*(2*q2_x.^2./(q2.^3) - q2_xx./(q2.^2));
s2_xxx = p2*(-6*q2_x.^3./(q2.^4) + 6*q2_x.*q2_xx./(q2.^3) - q2_xxx./(q2.^2));

DS = s2 - s1;
DS_x = s2_x - s1_x;
DS_xx = s2_xx - s1_xx;
DS_xxx = s2_xxx - s1_xxx;

%% ----------------- Time Of Flight Solution --------------------------- %%
% if sim.TOF_imposed_flag == 1
%     %- Interpolationg function coefficient a - fsolve + Cavalieri-Simpson to
%     %  solve the integral
%     coeff_a0  = 0;
%     fun = @(a) find_a(a,x,psi,TOF,sim,s1,s1_x,s1_xx,s2,s2_x,s2_xx,delta1,d1_x,delta1_xx,delta2,delta2_x,delta2_xx);
% 
%     %options=optimoptions('fsolve', 'TolFun', 1e-8, 'TolX', 1e-8,'Display','off');
%     options=optimoptions('fsolve','Display','none');
%     [coef_a,fsolve_err] = fsolve(fun,coeff_a0,options); %% alternative: fzero 
%     % a  = find_a_prof(x,psi,TOF,sim,s1,s1_x,s1_xx,s2,s2_x,s2_xx,delta1,delta1_x,delta1_xx,delta2,delta2_x,delta2_xx); % prof
% else
%     coef_a = 0;
% end
coef_a = 0;
coef_b = -20 - 4*coef_a;
coef_c = 70 + 6*coef_a;
coef_d = - 84 - 4*coef_a;
coef_e = 35 + coef_a;

%- Interpolating function xi(x,a) and derivatives
% Initial and final position and velocity shall be constrained, therefore
% the function shall be continuous (between 0 and 1).
% If you want TOF free (not imposed) : a = 0.     
X     = coef_a*x.^8 + coef_b*x.^7 + coef_c*x.^6 + coef_d*x.^5 + coef_e*x.^4;
X_x   = 8*coef_a*x.^7 + 7*coef_b*x.^6 + 6*coef_c*x.^5 + 5*coef_d*x.^4 + 4*coef_e*x.^3;
X_xx  = 56*coef_a*x.^6 + 42*coef_b*x.^5 + 30*coef_c*x.^4 + 20*coef_d*x.^3 + 12*coef_e*x.^2;
X_xxx = 336*coef_a*x.^5 + 210*coef_b*x.^4 + 120*coef_c*x.^3 + 60*coef_d*x.^2 + 24*coef_e*x;

%- Sun distance s(x) and declination delta(x) 
s     = DS.*X + s1;
s_x   = DS_x.*X + DS.*X_x + s1_x;
s_xx  = DS_xx.*X + 2*DS_x.*X_x + DS.*X_xx + s1_xx;
s_xxx = DS_xxx.*X + X_xxx.*DS + 3*DS_xx.*X_x + 3*DS_x.*X_xx + s1_xxx;

d = Dd.*X + d1;
d_x   = Dd_x.*X + Dd.*X_x + d1_x;
d_x_2 = d_x.^2;
d_x_3 = d_x.^3;
d_xx  = Dd_xx.*X + Dd.*X_xx + 2*Dd_x.*X_x + d1_xx;
d_xxx = Dd_xxx.*X + X_xxx.*Dd +3*Dd_xx.*X_x + 3*Dd_x.*X_xx + d1_xxx;

A1 = s_x;
B1 = s.*d_x;
A2 = s_xx - s.*d_x_2;
B2 = 2*s_x.*d_x + s.*d_xx;

%- Transformation between spherical and cylindrical coordinates
r = s.*cos(d);
r_x   = A1.*cos(d) - B1.*sin(d);
r_xx  = A2.*cos(d) - B2.*sin(d);
r_xxx = - (3*(s_xx.*d_x + s_x.*d_xx) + s.*(d_xxx - d_x_3)).*sin(d) + ...
          (s_xxx - 3*d_x.*(s_x.*d_x + s.*d_xx)).*cos(d);

z = s.*sin(d);
z_x   = A1.*sin(d) + B1.*cos(d);
z_xx  = B2.*cos(d) + A2.*sin(d);

%% Parametrized EOM
tan_gamma = r_x./(PSY*r);
s_3 = s.^3;

% Quantities to compute the second derivative of x wrt time
Nu = sim.mu*r; 
De = s_3.*(PSY_2.*r - r_xx + 2*PSY.*tan_gamma.*r_x); 

Nu_x = sim.mu*r_x;
De_x = 3*s_x./s.*De + s_3.*(r_x*PSY_2 - r_xxx + (2*r.*r_x.*r_xx - r_x.^3)./(r.^2));

% Square of derivative of x wrt time : x_t_2
x_t_2 = Nu./De; % \dot{x}^2
x_t   = sqrt(x_t_2);
x_tt = 0.5*(Nu_x - x_t_2.*De_x)./De;

%% TOF
% step zero of newton method
dx = x(2) - x(1);
three_odd_n_solm2 = 3:2:sim.n_sol-2;
two_even_n_solm1 = 2:2:sim.n_sol-1;

d_time_k1 = 1./x_t(1);
d_time_knumel_k = 1./x_t(sim.n_sol);
d_time_m = 1./x_t(two_even_n_solm1);
d_time_ktwo_end_m1 = 1./x_t(three_odd_n_solm2);
TOFc = (d_time_k1 + d_time_knumel_k)/6 + sum(d_time_ktwo_end_m1)/3;
TOFc = dx*2*(TOFc + 2/3*sum(d_time_m));
error = (TOF-TOFc)/TOF;

coef_a_old = coef_a;
coef_a = 1e-5;
count = 0;

while abs(error) > 1e-6 && count <= 100
    coef_b = -20 - 4*coef_a;
    coef_c = 70 + 6*coef_a;
    coef_d = - 84 - 4*coef_a;
    coef_e = 35 + coef_a;
    
    X     = coef_a*x.^8 + coef_b*x.^7 + coef_c*x.^6 + coef_d*x.^5 + coef_e*x.^4;
    X_x   = 8*coef_a*x.^7 + 7*coef_b*x.^6 + 6*coef_c*x.^5 + 5*coef_d*x.^4 + 4*coef_e*x.^3;
    X_xx  = 56*coef_a*x.^6 + 42*coef_b*x.^5 + 30*coef_c*x.^4 + 20*coef_d*x.^3 + 12*coef_e*x.^2;
    X_xxx = 336*coef_a*x.^5 + 210*coef_b*x.^4 + 120*coef_c*x.^3 + 60*coef_d*x.^2 + 24*coef_e*x;

    s     = DS.*X + s1;
    s_x   = DS_x.*X + DS.*X_x + s1_x;
    s_xx  = DS_xx.*X + 2*DS_x.*X_x + DS.*X_xx + s1_xx;
    s_xxx = DS_xxx.*X + X_xxx.*DS + 3*DS_xx.*X_x + 3*DS_x.*X_xx + s1_xxx;

    d = Dd.*X + d1;
    d_x   = Dd_x.*X + Dd.*X_x + d1_x;
    d_x_2 = d_x.^2;
    d_xx  = Dd_xx.*X + Dd.*X_xx + 2*Dd_x.*X_x + d1_xx;
    d_xxx = Dd_xxx.*X + X_xxx.*Dd +3*Dd_xx.*X_x + 3*Dd_x.*X_xx + d1_xxx;

    A1 = s_x;
    B1 = s.*d_x;
    A2 = s_xx - s.*d_x_2;
    B2 = 2*s_x.*d_x + s.*d_xx;

    r = s.*cos(d);
    r_x   = A1.*cos(d) - B1.*sin(d);
    r_xx  = A2.*cos(d) - B2.*sin(d);
    r_xxx = - (3*(s_xx.*d_x + s_x.*d_xx) + s.*(d_xxx - d_x_3)).*sin(d) + ...
          (s_xxx - 3*d_x.*(s_x.*d_x + s.*d_xx)).*cos(d);

    z = s.*sin(d);
    z_x   = A1.*sin(d) + B1.*cos(d);
    z_xx  = B2.*cos(d) + A2.*sin(d);
    
    tan_gamma = r_x./(PSY*r);
    s_3 = s.^3;

    Nu = sim.mu*r; 
    De = s_3.*(PSY_2.*r - r_xx + 2*PSY.*tan_gamma.*r_x); 

    Nu_x = sim.mu*r_x;
    De_x = 3*s_x./s.*De + s_3.*(r_x*PSY_2 - r_xxx + (2*r.*r_x.*r_xx - r_x.^3)./(r.^2));

    x_t_2 = Nu./De; % \dot{x}^2
    x_t   = sqrt(x_t_2);
    x_tt = 0.5*(Nu_x - x_t_2.*De_x)./De;
    
    % TOF
    d_time_k1 = 1./x_t(1);
    d_time_knumel_k = 1./x_t(sim.n_sol);
    d_time_m = 1./x_t(two_even_n_solm1);
    d_time_ktwo_end_m1 = 1./x_t(three_odd_n_solm2);
    % integral of time derivative
    TOFc = (d_time_k1 + d_time_knumel_k)/6 + sum(d_time_ktwo_end_m1)/3;
    TOFc = dx*2*(TOFc + 2/3*sum(d_time_m));
    error_old = error;
    error = (TOF-TOFc)/TOF;
    
    % newton formulas method
    d_error = (error - error_old)/(coef_a - coef_a_old);
    
    coef_a_old = coef_a;
    coef_a = coef_a_old - error/d_error;
    
    count = count + 1;
end



% ----------------------------------------------------------------------- %
if isreal(coef_a)
%- Thrust angle gamma -> It is imposed equal to the flight path angle in
%  order to have the in-plane thrust imposed as tangential only. In this way
%  the x time derivative can be analytically computed removing the dependency
%  from thrust per unit mass. 
%gamma = atan(r_x./(r*psi));
gamma = atan(tan_gamma); % || atan

%- In plane thrust per unit mass
Tin2m  = 1./cos(gamma) .* (2*PSY*r_x.*x_t_2 + r.*PSY.*x_tt);

%- Out of plane thrust per unit mass
Tout2m = z_xx.*x_t_2 + z_x.*x_tt + sim.mu./(s.^3).*z;

%- Thrust (absolute value) per unit mass
T2m = sqrt(Tin2m.^2 + Tout2m.^2);

%- Mass time history from Tsiolkovsky equation: m_t
   % devo scriverle in funzione di theta o x (||) ???

theta = Px;
dtheta = Px(2)- Px(1);
dtheta2 = dtheta/2;
dtheta12 = dtheta/12;
dtheta24 = dtheta/24;

theta_t = PSY*x_t;  

   if sim.direction == -1   %backward case (from dry to wet)
       m = theta; % o x? %% giusto?
       m(n_sol) = M;
       K = - T2m./(Isp*sim.g0)./theta_t; % dividi per theta_t per dimensionalità
       m(n_sol-1) = m(n_sol) -   dtheta2*(K(n_sol)*m(sim.n_sol)...
                        +K(n_sol-1)*m(n_sol-1)); % (?) Alternative: FE
       m(n_sol-2) = m(n_sol-1) - dtheta2*(K(n_sol-1)*m(n_sol-1)...
                        +K(n_sol-2)*m(n_sol-2));
        
       for i = n_sol-2:-1:2 % predictor corrector AB3AM4 
           m(i-1) = m(i) - dtheta12*(23*K(i)*m(i) - 16*K(i+1)*m(i+1) + 5*K(i+2)*m(i+2) ) ;
           m(i-1) = m(i) - dtheta24*(9*K(i-1)*m(i-1)+19*K(i)*m(i) - 5*K(i+1)*m(i+1) + K(i+2)*m(i+2) ) ;
       end
       
   else %forward case (from wet to dry)
        
        m = theta; %% o x? giusto?
        m(1) = M;
        K = - T2m./(Isp*sim.g0)./theta_t;
        m(2) = m(1) + dtheta2*( K(1)*m(1)+K(2)*m(2) );
        m(3) = m(2) + dtheta2*( K(2)*m(2)+K(3)*m(3) );
        
        for i = 3:n_sol-1 % dal quarto al penultimo punto con predictor corrector AB3AM4 expl
            m(i+1) = m(i) + dtheta12*( 23*K(i)*m(i) - 16*K(i-1)*m(i-1) + 5*K(i-2)*m(i-2) );
            m(i+1) = m(i) + dtheta24*( 9*K(i+1)*m(i+1)+19*K(i)*m(i) - 5*K(i-1)*m(i-1) + ...
                K(i-2)*m(i-2) );
        end
        
    end

   
% time vector
d_time =1./x_t;
dx = x(2) - x(1);

t    = x;
t(1) = 0; % Initialization

for i = 2:n_sol-1
    t(i) = t(i-1) + dx/6 * (d_time(i-1) + 4*d_time(i) + d_time(i+1) ) ; %% RK2 (?)
end
t(n_sol) = t(n_sol-1) + d_time(n_sol-1)*dx ;

% Thrust
Tin = Tin2m.*m * 1000* DU/TU^2;
Tout = Tout2m.*m * 1000* DU/TU^2;

    % Output
    output.m       = m ;
    output.t       = t ;
    output.Thrust  = [Tin,gamma,Tout];
    output.r       = r;
    output.theta   = theta;
    output.z       = z;
    output.a       = coef_a;
    output.Href    = h_ref_vers;

else
    % Output in case a is complex -> no solution with that tof
    output.m       = nan ;
    output.t       = nan ;
    output.Thrust  = nan ;
    output.r       = nan ;
    output.theta   = nan ;
    output.z       = nan ;
    output.a       = nan ;
    output.Href    = nan ;
end
end