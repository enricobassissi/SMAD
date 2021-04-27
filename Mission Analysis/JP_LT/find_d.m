function [ error_TOF ] = find_d( plan_d , plan_a , plan_b , plan_c  , theta_f , rf , tan_gam_f , rate_theta_f , TOF,l,l2,l3,l4,l5,l6,n_sol)

% calculation of e f g
AAA = [  30*theta_f^2 , -10*theta_f^3 ,    theta_f^4  ;
    -48*theta_f   ,  18*theta_f^2 , -2*theta_f^3  ;
    20           ,  -8*theta_f   ,    theta_f^2 ];

bbb = [ 1/rf - (plan_a + plan_b*theta_f + plan_c*theta_f^2 + plan_d*theta_f^3 );
    -(tan_gam_f)/rf - (plan_b + 2*plan_c*theta_f + 3*plan_d*theta_f^2);
    1/((rf^4)*(rate_theta_f^2)) - (1/rf + 2*plan_c + 6*plan_d*theta_f)];

e = 0.5/theta_f^6*AAA(1,:)*bbb;
f = 0.5/theta_f^6*AAA(2,:)*bbb;
g = 0.5/theta_f^6*AAA(3,:)*bbb;

lm = 0.5*(l(1:end-1) + l(2:end));

lm2 = lm.*lm;
lm3 = lm2.*lm;
lm4 = lm3.*lm;
lm5 = lm4.*lm;
lm6 = lm5.*lm;
% calculation r at each theta and theta_m
r  = 1./(plan_a + plan_b*l  + plan_c*l2  + plan_d*l3  + e*l4  + f*l5  + g*l6) ;
rm = 1./(plan_a + plan_b*lm + plan_c*lm2 + plan_d*lm3 + e*lm4 + f*lm5 + g*lm6) ;

% calculation of dt at each theta
dTOF  = sqrt( abs( r.^4  .* (1./r  + 2*plan_c + 6*plan_d*l  + 12*e*l2  + 20*f*l3  + 30*g*l4))  ) ;
dTOFm = sqrt( abs( rm.^4 .* (1./rm + 2*plan_c + 6*plan_d*lm + 12*e*lm2 + 20*f*lm3 + 30*g*lm4)) ) ;

% step size
h = l(2)-l(1);

% Cavalieri method
I = 0;
for i = 2:n_sol
    I = I + dTOF(i-1) + 4*dTOFm(i-1) + dTOF(i);
end

I = I*h/6;

% computation of TOF error (residual)
error_TOF =TOF-I;


end

