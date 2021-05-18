%% validation of delta az, el for cone of vinf
clear
clc

r = [1,1.2,0.1]; 
v = [1.1, 2.2, 0.01]; 
vv = v./norm(v); 
% [azimuthr,elevationr,radiusr] = cart2sph(r(1),r(2),r(3)); 
% translation on the point of application
[azimuthv,elevationv,radiusv] = cart2sph(v(1),v(2),v(3)); 

azimuthv2 = azimuthv + pi/30;
elevationv2 = elevationv + pi/30;
[v2x,v2y,v2z] = sph2cart(azimuthv2,elevationv2,radiusv);
v2 = [v2x,v2y,v2z];

azimuthv3 = azimuthv + pi/9;
elevationv3 = elevationv + pi/9;
[v3x,v3y,v3z] = sph2cart(azimuthv3,elevationv3,radiusv);
v3 = [v3x,v3y,v3z];

azimuthv4 = azimuthv + pi/3;
elevationv4 = elevationv + pi/3;
[v4x,v4y,v4z] = sph2cart(azimuthv4,elevationv4,radiusv);
v4 = [v4x,v4y,v4z];

figure()
plot3(r(1),r(2),r(3),'*')
hold on
plot3(0,0,0,'m*')
quiver3(r(1),r(2),r(3),v(1),v(2),v(3),'DisplayName','v')
quiver3(r(1),r(2),r(3),v2(1),v2(2),v2(3),'DisplayName','v2')
quiver3(r(1),r(2),r(3),v3(1),v3(2),v3(3),'DisplayName','v3')
quiver3(r(1),r(2),r(3),v4(1),v4(2),v4(3),'DisplayName','v4')
legend('show')
axis equal

%% aim cone
clear
clc

r = [1,1.2,0.1]; 
v = [1.1, 2.2, 0.01]; 

% [azimuthr,elevationr,radiusr] = cart2sph(r(1),r(2),r(3)); 
% translation on the point of application
% [azimuthv,elevationv,radiusv] = cart2sph(r(1)+v(1),r(2)+v(2),r(3)+v(3)); 
[azimuthv,elevationv,radiusv] = cart2sph(v(1),v(2),v(3)); 

figure()
plot3(r(1),r(2),r(3),'*')
hold on
plot3(0,0,0,'m*')
quiver3(r(1),r(2),r(3),v(1),v(2),v(3),'--','DisplayName','v')
lim_sotto = -deg2rad(20);
lim_sopra = deg2rad(20);
step = deg2rad(10);
for i = lim_sotto:step:lim_sopra
    for j = lim_sotto:step:lim_sopra
        azimuthv2 = azimuthv + i;
        elevationv2 = elevationv + j;
        [v2x,v2y,v2z] = sph2cart(azimuthv2,elevationv2,radiusv);
        v2 = [v2x,v2y,v2z];

        quiver3(r(1),r(2),r(3),v2(1),v2(2),v2(3),'DisplayName',num2str([i,j]))
    end
end
legend('show')
axis equal