function f = dynamic(x,u)
mu=1; %adim
% [r(1),r(2),r(3)]=sph2cart(x(1),x(2),x(3));
% g=-mu*r./(norm(r)^3);
% [v(1),v(2),v(3)]=Vspher2car(x(1),x(2),x(3),x(4),x(5),x(6));
% [g_my(1),g_my(2),g_my(3)]=Acart2sph(r(1),r(2),r(3),v(1),v(2),v(3),g(1),g(2),g(3));
r_=x(2);
g=-mu/(r_^2);

f=[x(3)/r_; x(4); -(x(4)*x(3)/r_)+u(1); x(3)^2/r_+g+u(2)];
end