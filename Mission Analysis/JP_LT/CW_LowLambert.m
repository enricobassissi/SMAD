function [ output] = CW_LowLambert( RI , RF , VI , VF , N_rev , TOF ,M ,hp , kp , PS ,sim )

% Measurament units definition
DU = sim.DU ;                    % Distance unit lenght[km]
TU = sim.TU;                     % Time unit duration [s]
n_sol = sim.n_sol;               % number of computational nodes

% magnitude of the Sun distance
ri = sqrt(RI(1)*RI(1) + RI(2)*RI(2) + RI(3)*RI(3));
rf = sqrt(RF(1)*RF(1) + RF(2)*RF(2) + RF(3)*RF(3));

% initial and final position vectors
RI_vers = RI/ri;
RF_vers = RF/rf;

% scalar and dot product between initial and final versors
RIvcRFv = [ RI_vers(2)*RF_vers(3) - RI_vers(3)*RF_vers(2);
    RI_vers(3)*RF_vers(1) - RI_vers(1)*RF_vers(3) ;
    RI_vers(1)*RF_vers(2) - RI_vers(2)*RF_vers(1) ];

rivDrfv = RI_vers(1)*RF_vers(1) + RI_vers(2)*RF_vers(2) + RI_vers(3)*RF_vers(3);

% transfer plane definition (z = 0);
norm_RCRRv = sqrt(RIvcRFv(1)*RIvcRFv(1) + RIvcRFv(2)*RIvcRFv(2) + RIvcRFv(3)*RIvcRFv(3));
RCRRv = RIvcRFv/norm_RCRRv;

% the transfer must be counterclockwise
if RCRRv(3) < 0
    RCRRv = -RCRRv;
end

% determination of transfer angle
if norm_RCRRv < 1e-10 % singular case
    
    RCRRv = [ RI_vers(2)*VI(3) - RI_vers(3)*VI(2);
        RI_vers(3)*VI(1) - RI_vers(1)*VI(3);
        RI_vers(1)*VI(2) - RI_vers(2)*VI(1)];
    
    RCRRv = RCRRv/sqrt(RCRRv(1)*RCRRv(1) + RCRRv(2)*RCRRv(2) + RCRRv(3)*RCRRv(3));
    
    if rivDrfv > 0
        psy = 0;
    else
        psy = pi;
    end
else
    if RIvcRFv(3) > 0 % general case
        psy = acos(rivDrfv);
    else
        psy = 2*pi - acos(rivDrfv);
    end
end

% compute final anomaly adding numer of revolution
L = psy + 2*N_rev*pi ;

L2 = L*L;
L3 = L2*L;
L4 = L3*L;
L5 = L4*L;
L6 = L5*L;

% initial and final velocity decomposition
v_perp_i = VI(1)*RCRRv(1) + VI(2)*RCRRv(2) + VI(3)*RCRRv(3);
V_PLAN_I = VI-v_perp_i*RCRRv;


v_perp_f = VF(1)*RCRRv(1) + VF(2)*RCRRv(2) + VF(3)*RCRRv(3);
V_PLAN_F = VF-v_perp_f*RCRRv;


% -------------- PLANAR MOTION --------------
% in plane initial velocity decomposition
VI_rad = RI_vers*(V_PLAN_I(1)*RI_vers(1) + V_PLAN_I(2)*RI_vers(2) + V_PLAN_I(3)*RI_vers(3));
VI_tra = V_PLAN_I-VI_rad;
norm_VI_tra = sqrt(VI_tra(1)*VI_tra(1) + VI_tra(2)*VI_tra(2) + VI_tra(3)*VI_tra(3));

% in plane final velocity decomposition
VF_rad =  RI_vers*(V_PLAN_F(1)*RF_vers(1) + V_PLAN_F(2)*RF_vers(2) + V_PLAN_F(3)*RF_vers(3));
VF_tra = V_PLAN_F-VF_rad;
norm_VF_tra = sqrt(VF_tra(1)*VF_tra(1) + VF_tra(2)*VF_tra(2) + VF_tra(3)*VF_tra(3));

% initial and final flight path angles
tan_gam_i = (V_PLAN_I(1)*RI_vers(1) + V_PLAN_I(2)*RI_vers(2) + V_PLAN_I(3)*RI_vers(3))/norm_VI_tra;
tan_gam_f = (V_PLAN_F(1)*RF_vers(1) + V_PLAN_F(2)*RF_vers(2) + V_PLAN_F(3)*RF_vers(3))/norm_VF_tra;

% initial and final angular velocities
rate_theta_i = norm_VI_tra/ri;
rate_theta_f = norm_VF_tra/rf;

% computation of Conway inverse polinomium coefficients a,b,c
plan_a = 1/ri;
plan_b = -tan_gam_i / ri ;
plan_c = (1/(2*ri)) * ((1/( (ri^3) * rate_theta_i^2 )) - 1) ;

% generating theta vector for the solution
l  = L*sim.x;
l2 = l.*l;
l3 = l2.*l;
l4 = l3.*l;
l5 = l4.*l;
l6 = l5.*l;

% solving for d
fun = @(plan_d) find_d( plan_d , plan_a , plan_b , plan_c , L , rf , tan_gam_f , rate_theta_f , TOF ,l,l2,l3,l4,l5,l6,sim.n_sol);
plan_d = 0;
% opt = optimset('Display','off');
% plan_d = fzero(fun,plan_d,opt);
plan_d = JP_secant_solver( fun , plan_d , 1e-8 , 1e-6 , 100);

% calculation of e f g
AAA = [  30*L2 , -10*L3 ,    L4  ;
    -48*L   ,  18*L2 , -2*L3  ;
    20      ,  -8*L   ,    L2 ];

bbb = [ 1/rf - (plan_a + plan_b*L + plan_c*L2 + plan_d*L3 );
    -(tan_gam_f)/rf - (plan_b + 2*plan_c*L + 3*plan_d*L2);
    1/((rf^4)*(rate_theta_f^2)) - (1/rf + 2*plan_c + 6*plan_d*L)];

plan_e = 0.5/L6*AAA(1,:)*bbb;
plan_f = 0.5/L6*AAA(2,:)*bbb;
plan_g = 0.5/L6*AAA(3,:)*bbb;


% radius at each theta
r = 1./(plan_a + plan_b*l + plan_c*l2 + plan_d*l3 + plan_e*l4 + plan_f*l5 + plan_g*l6) ;
r2 = r.*r;
r3 = r.*r2;
r4 = r.*r3;

% angoular velocity
l_d = sqrt((1./r4)./(1./r + 2*plan_c + 6*plan_d*l + 12*plan_e*l2 + 20*plan_f*l3 + 30*plan_g*l4));

if isreal(l_d)
    
    % fligth path angle at each theta
    tan_gam = -r.*(plan_b + 2*plan_c*l + 3*plan_d*l2 + 4*plan_e*l3 + 5*plan_f*l4 + 6*plan_g*l5);
    gam = atan(tan_gam);
    
    %rrt
    % rrt
    l_dd = -1./(2.*r4).*( (4*tan_gam)./((1./r)+2*plan_c + 6*plan_d*l + 12*plan_e*l2 + 20*plan_f*l3 + 30*plan_g*l4) + (6*plan_d+24*plan_e*l+60*plan_f*l2 + 120*plan_g*l3 - tan_gam./r)./((1./r)+2*plan_c + 6*plan_d*l + 12*plan_e*l2 + 20*plan_f*l3 + 30*plan_g*l4).^2 );
    
    % radial velocity
    v_r = -r2.*l_d.*(plan_b + 2*plan_c*l + 3*plan_d*l2 + 4*plan_e*l3 + 5*plan_f*l4 + 6*plan_g*l5);
    
    % adimensional and dimensional acceleration
    plane_acc = -sim.mu./(2*r3.*cos(gam)).*(6*plan_d+24*plan_e*l+60*plan_f*l2+120*plan_g*l3-tan_gam./r)./(1./r + 2*plan_c + 6*plan_d*l + 12*plan_e*l2 + 20*plan_f*l3 + 30*plan_g*l4).^2;
    
    % time vector
    d_time =1./l_d;
    dl = l(2)-l(1);
    t = l;
    
    t(1) = 0;
    for i = 2:n_sol-1
        t(i) = t(i-1) + dl/6 * (d_time(i-1) + 4*d_time(i) + d_time(i+1) ) ;
    end
    t(end) = t(end-1) + d_time(end-1)*dl ;
    
    
    % OUT OF PLANE MOTION
    
    switch sim.out_shape
        
        case 0 % planar motion only
            
            z = 0*l;
            s = r;
            v_z = z;
            z_dd = z;
            
            
        case 1 % J shape
            
            a11 = L*exp(-kp*L);
            a12 = L*exp(-hp*L);
            a13 = L ;
            a21 = l_d(1);
            a22 = (1-hp*L)*exp(-hp*L)*l_d(1);
            a23 = l_d(1);
            a31 = (1-kp*L)*exp(-kp*L)*l_d(end);
            a32 = l_d(end);
            a33 = l_d(end);
            
            detA = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a31*a22*a13 - a32*a23*a11 - a21*a12*a33 ;
            
            invA = 1/detA * [a22*a33-a23*a32 , a13*a32-a12*a33 , a12*a23-a13*a22 ;
                a23*a31-a21*a33 , a11*a33-a13*a31 , a13*a21-a11*a23 ;
                a21*a32-a22*a31 , a12*a31-a11*a32 , a11*a22-a12*a21 ];
            
            out = invA * [0;v_perp_i;v_perp_f];
            
            out_a = out(1);
            out_b = out(2);
            out_c = out(3);
            out_d = out_b*L*exp(-hp*L);
            
            z = out_a*l.*exp(-kp*l) + out_b*(l-L).*exp(hp*(l-L)) + out_c*l + out_d;
            
            s = sqrt(r2+z.^2);
            
            v_z =  l_d.*( out_a*exp(-kp*l) - out_a*kp*l.*exp(-kp*l) + out_b*exp(hp*(l-L)) + out_b*hp*(l-L).*exp(hp*(l-L))   + out_c );
            
            z_dd = l_dd.*( out_a*exp(-kp*l) - out_a*kp*l.*exp(-kp*l) + out_b*exp(hp*(l-L)) + out_b*hp*(l-L).*exp(hp*(l-L))   + out_c ) + l_d.^2.*( out_a*kp*exp(-kp*l).*(kp*l-2) + out_b*hp*exp(hp*(l-L)).*(hp*(l-L)+2) );
            
        case 2 % Conway-Wall shape
            
            q = hp;
            
            b_z = v_perp_i/rate_theta_i;
            
            temp = 1/L^q * [q*L , -L^2 ; -(q-1) , L] * [-b_z*L ; v_perp_f/rate_theta_f - b_z];
            c_z = temp(1);
            d_z = temp(2);
            
            z    = b_z*l + c_z*l.^(q-1) + d_z*l.^q;
            v_z  = l_d.*(b_z + (q-1)*c_z.*l.^(q-2) + q*d_z*l.^(q-1) );
            z_dd = l_dd.*(b_z + (q-1)*c_z.*l.^(q-2) + q*d_z*l.^(q-1) ) + l_d.^2.*((q-1)*(q-2)*c_z.*l.^(q-3) + q*(q-1)*d_z*l.^(q-2) );
            s = sqrt(r2+z.^2);
    end
    
    out_acc = z_dd + sim.mu*z./s.^3;
    
    %thrust and mass profile
    a_vect = sqrt(out_acc.^2 + plane_acc.^2);
    
    dl = l(2)-l(1);
    dl12 = dl/12;
    dl24 = dl/24;
    
    if sim.direction == -1
        
        m = l;
        m(sim.n_sol) = M;
        K = -a_vect/(PS.Is*sim.g0)./l_d;
        m(sim.n_sol-1) = m(sim.n_sol) -dl/2*(K(sim.n_sol)*m(sim.n_sol)+K(sim.n_sol-1)*m(sim.n_sol-1));
        m(sim.n_sol-2) = m(sim.n_sol-1) -dl/2*(K(sim.n_sol-1)*m(sim.n_sol-1)+K(sim.n_sol-2)*m(sim.n_sol-2));
        
        for i = sim.n_sol-2:-1:2 % dal quarto al penultimo punto con predictor corrector AB3AM4 expl
            m(i-1) = m(i) - dl12*(23*K(i)*m(i) - 16*K(i+1)*m(i+1) + 5*K(i+2)*m(i+2) ) ;
            m(i-1) = m(i) - dl24*(9*K(i-1)*m(i-1)+19*K(i)*m(i) - 5*K(i+1)*m(i+1) + K(i+2)*m(i+2) ) ;
        end
        
    else
        
        m = k;
        m(1) = M;
        K = -a_vect/(PS.Is*sim.g0)./l_d;
        m(2) = m(1) +dl/2*(K(1)*m(1)+K(2)*m(2));
        m(3) = m(2) +dl/2*(K(2)*m(2)+K(3)*m(3));
        
        for i = 3:sim.n_sol-1 % dal quarto al penultimo punto con predictor corrector AB3AM4 expl
            m(i+1) = m(i) + dl12*(23*K(i)*m(i) - 16*K(i-1)*m(i-1) + 5*K(i-2)*m(i-2) ) ;
            m(i+1) = m(i) + dl24*(9*K(i+1)*m(i+1)+19*K(i)*m(i) - 5*K(i-1)*m(i-1) + K(i-2)*m(i-2) ) ;
        end
        
    end
    
    
    plane_T = plane_acc.*m * 1000* DU/TU^2;
    out_T = out_acc.*m * 1000* DU/TU^2;
    
    output.u       = [plane_T,gam,out_T] ;
    output.m       = m ;
    output.t       = t ;
    output.s       = s ;
    output.r       = r ;
    output.l       = l ;
    output.z       = z ;
    output.RI      = RI ;
    output.RF      = RF ;
    output.h_ref_v = RCRRv ;
    output.PSY     = L ;
    output.v_r     = v_r ;
    output.w       = l_d;
    output.v_z     = v_z ;
    
else
    
    output.u       = nan ;
    output.m       = nan ;
    output.t       = nan ;
    output.s       = nan ;
    output.r       = nan ;
    output.Px      = nan ;
    output.z       = nan ;
    output.RI      = nan ;
    output.RF      = nan ;
    output.h_ref_v = nan ;
    output.PSY     = nan ;
    output.x_t     = nan ;
    output.z_d     = nan ;
    
end

end
