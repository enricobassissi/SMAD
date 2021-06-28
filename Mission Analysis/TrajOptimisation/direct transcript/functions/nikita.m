%% ----- stuff for nikita
function PROP = nikita(SC, data, sim, colors)

N = data.n_int;
tol = 1.5e-3;
TC = SC.HS.T > tol;
% figure('Name','T')
% plot(SC.timead*sim.TU/86400, SC.HS.T,'Color', colors(2,:));
% xlabel('t [d]'); ylabel('T [N]');
% title('Thrust');
% hold on
% yline(tol)

% ----- how many duty cycles
if TC(1) == 1
	DutyCycle = 1;
else
   DutyCycle = 0;
end
for i=2:length(TC)
    if TC(i)-TC(i-1)>0
        DutyCycle = DutyCycle + 1;
    end
end
PROP.DutyCycle = DutyCycle;

% ---- firing time
StepTime = SC.timead(2) - SC.timead(1);
PROP.ThrustTime = sum(StepTime*sim.TU/86400*TC);

% ---- total impulse
% ---- integral of thrust on time
PROP.TotalImpulse = trapz(SC.timead*sim.TU,SC.HS.T); % Ns

% --- distance from the sun
PROP.DistanceSC_Sun_magnitude = SC.HS.r;

% ---- mass after the leg, at beginning of coasting
PROP.mass_end = SC.HS.X(end,7); % last node, 7th param is the mass

% --- POWER modulation over distance, time (degradation), angle wrt sun
    switch data.ThrustModulationFlag
        case 0
            eta = 1;
        case 1
            % distance
            eta_space = solar_irradiance_cooler(SC.HS.r, 1, 1); % r is in AU
            eta = eta_space/max(eta_space);
        case 2
            % distance
            eta_space = solar_irradiance_cooler(SC.HS.r, 1, 1); % r is in AU
            
            % degradation
            % REFERENCES:
            D = 0.01674;
            time_in_years = SC.timead*sim.TU/(86400*365);
            eta_degradation = (1-D).^time_in_years;

            %total efficency
            eta = eta_space.*eta_degradation;
            eta = eta/max(eta);
        case 3
            % distance
            eta_space = solar_irradiance_cooler(SC.HS.r, 1, 1); % r is in AU
            
            % degradation
            % REFERENCES:
            D = 0.01674;
            time_in_years = SC.timead*sim.TU/(86400*365);
            eta_degradation = (1-D).^time_in_years;

            % angle of panels respect to sun
            % --- local cartesian
            [xx,yy,zz] = pol2cart(SC.HS.X(:,2),SC.HS.X(:,1),SC.HS.X(:,3));
            in_plane_position = [xx,yy,zeros(N,1)];
            in_plane_position_norm = vecnorm(in_plane_position,2,2);
            in_plane_position_vers = in_plane_position./in_plane_position_norm;
        
            z_vers = [0,0,1].*ones(N,3);
            tg_vers = -cross(in_plane_position_vers,z_vers);
        
            % ----- cartesian 3D thrust
            % beta = elevation, alpha = thrust angle in plane
            % x -> T t , y = -T rad, z -> T outofplane
            [Tx, Ty, Tz] = sph2cart(SC.HS.alpha',SC.HS.beta',SC.HS.T');
            % cartesiano centrato movente con la spacecraft
            Ttg = Tx;
            Tradial = -Ty;
            Tout = Tz;
        
            T_cart_radial = Tradial.*in_plane_position_vers; 
            T_cart_tg = Ttg.*tg_vers;
            T_cart_out = Tout.*z_vers;
        
            T_cart_local = T_cart_radial+T_cart_tg+T_cart_out;
        
            % check -- if low it's the cross product
        	% min(T' - vecnorm(T_cart_local,2,2))
        	% max(T' - vecnorm(T_cart_local,2,2))
        
            % -- T in cartesiano
            TT = rotate_local2ecplitic(SC.r_in,T_cart_local,N,SC.Href);
        
            r_transf_orbit  = [SC.HS.X(:,1).*cos(SC.HS.X(:,2)), SC.HS.X(:,1).*sin(SC.HS.X(:,2)), SC.HS.X(:,3)];
            RR = rotate_local2ecplitic(SC.r_in,r_transf_orbit,N,SC.Href);
        
        	% RR,TT, must be Nx3, cartesian position and thrust
            a_sun = acos(dot(RR,TT,2)./(vecnorm(TT,2,2).*SC.HS.X(:,1))) -pi/2; 
            eta_angle = cos(a_sun);
            eta_angle(isnan(eta_angle)) = cosd(4);

            %total efficency
            eta = eta_space.*eta_angle.*eta_degradation;
            eta = eta/max(eta);
    end
    
    % -- POWER modulation and checks
    P_modulated = data.Pmax*eta;
    
    % the thruster we selected
    [~, T_max_modulated, ~] = qinetiq_t5_future(P_modulated); % column vector
    
figure('Name','T')
plot(SC.timead, SC.HS.T,'Color', colors(1,:),'DisplayName','T requested')
hold on; grid on;
plot(SC.timead, T_max_modulated,'--','Color', colors(2,:),'DisplayName','T available')
legend('show')

end