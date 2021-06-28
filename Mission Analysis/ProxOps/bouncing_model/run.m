clear all
time_sim=10100;
w=50*2*pi/(24*3600);
m=1;
e_coeff=0.3;
r_bouncer=0.1;
I=eye(3).*(0.4*m*r_bouncer^2); %inertia matrix of the bouncer
f=0.04;

out=sim('bouncing_int.slx');
figure()
plot3(out.x_ast.Data(:),out.y_ast.Data(:),out.z_ast.Data(:));
hold on
[x,y,z] = sphere;
x = x*3.5;
y = y*3.5;
z = z*3.5;
surf(x,y,z)
axis equal
% figure()
% plot(out.x_ast.Data(:),out.y_ast.Data(:))
% hold on
% plot(3.5*cos(linspace(0,2*pi)),3.5*sin(linspace(0,2*pi)))
% figure()
% comet3(out.x_ast.Data(:),out.y_ast.Data(:),out.z_ast.Data(:));
