function [] = plot_orbit(z, parout, data)

%{
 Function that returns the following plots: (1) 3D orbit evolution, (2) Ground track 

 INPUT:  1. z: output of numerical integration (states evolution in time: [a e i OM om th xa va Vout Vout_int I xv vv]
         2. parout: additional output of numerical integration [Th D A_valve lon lat height r v at an ah]
         3. data: Characteristic data of GOCE 

 FUNCTIONS REQUIRED: plot_groundTrack, plot_earth

 CONTRIBUTORS:  Bassissi Enrico
                Colombo Alessandro
                De Luca Maria Alessandra

%}

% Extraction of quantities of interest
lon = parout(:,4);        % [deg]
lat = parout(:,5);        % [deg]

% Recovering of position and velocity 
coe_plot = [z(:,1) z(:,2) z(:,3) z(:,4) wrapTo2Pi(z(:,5)) wrapTo2Pi(z(:,6))];  
y_plot = zeros(size(coe_plot));
for k = 1:length(coe_plot)
    [r_plot,v_plot] = sv_from_coe(coe_plot(k,:),data.earth.mu_E);
    y_plot(k,:) = [r_plot', v_plot'];
end

% Plot orbit evolution
figure('Name', '3D Orbit')
plot_earth;
hold on
plot3 (y_plot(:,1), y_plot(:,2), y_plot(:,3),'b');
title('Orbit Evolution')
grid minor; xlabel('km'); ylabel('km'); zlabel('km')
view([45 45])

% Plot Ground track 
figure('Name', 'Ground Track') 
x=[-180 180];
y = [-90 90];
C=imread('EarthTexture.jpg');
C=flip(C,1);
image(x,y,C);
title('Ground Track')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
axis xy
hold on
h1 = plot_groundTrack(lon,lat,'g'); 


end

