function [SAA,EVA,SCA,SolarConjunction] = aspect_angles(sol)

    % SAA: Sun Aspect Angle
    % EVA: Earth Visibility Angle
    % SCA: Solar Conjunction Angle
    
    % position and velocity are actually already wrt sun
    position_wrt_sun = sol.SpacecraftTrajectory(:,1:3);
    velocity_wrt_sun = sol.SpacecraftTrajectory(:,4:6);
    
%     % EXAMPLE to explain method
%     a = 1e-10  % start with a very small angle
%     the differences start to happen under 1e-7
%     this method is resistant also to high angles, pi,3/2pi ecc
%     u = 4*[1 0 0]  % arbitrary non-unit vector in X direction
%     v = 5*[cos(a) sin(a) 0]  % vector different from u by small angle
%     acos(dot(u,v)/(norm(u)*norm(v))) % acos formulation does not recover the small angle
%     ans =0
%     atan2(norm(cross(u,v)),dot(u,v)) % atan2 formulation does recover the small angle
%     ans = 1e-10
    
    % Initialise the vector
    SAA = zeros(size(position_wrt_sun,1),1);
    EVA = zeros(size(position_wrt_sun,1),1);
    SCA = zeros(size(position_wrt_sun,1),1);
    SolarConjunction = zeros(size(position_wrt_sun,1),1);
    for i=1:size(position_wrt_sun,1)
        pos_vers = position_wrt_sun(i,:)./norm(position_wrt_sun(i,:));
        vel_vers = velocity_wrt_sun(i,:)./norm(velocity_wrt_sun(i,:));
        % Sun Angle
        SAA(i,1) = atan2(norm(cross(vel_vers,pos_vers)),dot(vel_vers,pos_vers));
        
        % Direction of the Earth for Earth Comms
        mjd2000 = sol.SCtime(i)/(3600*24);
        [kep_EA,ksun] = uplanet(mjd2000, 3);
        [rEA, ~] = sv_from_coe(kep_EA,ksun); % km, km/s
        DistSpacecraftEarth = position_wrt_sun(i,:) - rEA';
        SC_EA_vers = DistSpacecraftEarth./norm(DistSpacecraftEarth);
%         rEA_norm = rEA./norm(rEA);
        % MAYBE
        EVA(i,1) = atan2(norm(cross(vel_vers,SC_EA_vers)),dot(vel_vers,SC_EA_vers));
        
        % Earth Covered by the sun
        % Solar Conjunction
        % Sun disturbance radius considered 10 times its radius
        SunRadius = astroConstants(3)*10; %km
%         DistSpacecraftEarth = position_wrt_sun(i,:) - rEA';
        % tol is the angle given by the Sun disk at a given distance
        tol = atan2(norm(cross(position_wrt_sun(i,:),(position_wrt_sun(i,:)-SunRadius))),...
                dot(position_wrt_sun(i,:),(position_wrt_sun(i,:)-SunRadius)));
%         SC_EA_norm = DistSpacecraftEarth./norm(DistSpacecraftEarth);
        
        SCA(i,1) = atan2(norm(cross(SC_EA_vers,pos_vers)),dot(SC_EA_vers,pos_vers));
        % when does the conjunction happens?
        % if the angle between sc-earth/sc-sun is less than tol
        % and the distance sc-earth > sc-sun, meaning that's behind
        if abs(SCA(i,1)) < tol && norm(DistSpacecraftEarth) > norm(position_wrt_sun(i,:))
            SolarConjunction(i,1) = 1;
        end

    end

end