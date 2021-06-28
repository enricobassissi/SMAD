function [c, ceq] = EoM(X, data,sim,DT)
    muS = data.muS;
    time = data.time;
    N = data.n_int;
    r_in = DT.r_in;
    Href = DT.Href;
        
    x = zeros(N,7);
    m = zeros(1, N);
    T = m; ualpha = m; ubeta = m;
    EoMpolarAD=@(t_hand, x_hand, T_hand, alpha_hand, beta_hand, muS_hand, data_hand) EoMpolar(t_hand, x_hand, T_hand, alpha_hand, beta_hand, muS_hand, data_hand,sim);
    for k=1:N
        %adimensional variables
        x(k,:)      = X((k-1)*10 + 1:(k-1)*10+7);
        % dimensional variables
        T(k)        = X((k-1)*10+8); 
        ualpha(k)   = X((k-1)*10+9);
        ubeta(k)    = X((k-1)*10+10); 
        
        m(k) = x(k,7);
    end
 
    %initial and final states
%     xi              = data.xi;
%     xf              = data.xf;
%     BCi             = (xi(1:6) - x(1,1:6))*sim.DU;
%     BCf             = (xf(1:6) - x(end,1:6))*sim.DU;
    
%     Xprop           = zeros(N, 7);
%     Xprop(1,:)      = x(1,:); %first value equal
    
%     %RK4 forward integration (in each interval of delta t)
%     for k = 1:N-1
% %         xx = Xprop(k,:);
%         xx = x(k,:);
% 
%         hhh = time(k+1) - time(k);
%         K1 = EoMpolarAD(time(k), xx, T(k), ualpha(k), ubeta(k), muS, data);
%         K2 = EoMpolarAD(time(k)+ hhh/2, xx + hhh/2*K1, 0.5*(T(k+1) + T(k)), 0.5*(ualpha(k+1) + ualpha(k)), 0.5*(ubeta(k)+ubeta(k+1)), muS, data);
%         K3 = EoMpolarAD(time(k)+ hhh/2, xx + hhh/2*K2, 0.5*(T(k+1) + T(k)), 0.5*(ualpha(k+1) + ualpha(k)), 0.5*(ubeta(k)+ubeta(k+1)),  muS, data);
%         K4 = EoMpolarAD(time(k)+ hhh, xx + hhh*K3, T(k+1), ualpha(k+1), ubeta(k+1), muS, data);    
% 
%         Xprop(k+1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
% %         Xprop(k+1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
%     end
%     
% %RK4 backward integration (in each interval of delta t)
%     Xprop(end,:)      = x(end,:); %first value equal
% 
%     for k = N:-1:2
%         xx = Xprop(k,:);
% %         xx = x(k,:);
% 
%         hhh = time(k-1) - time(k);
%         K1 = EoMpolarAD(time(k), xx, T(k), ualpha(k), ubeta(k), muS, data);
%         K2 = EoMpolarAD(time(k)+ hhh/2, xx + hhh/2*K1, 0.5*(T(k-1) + T(k)), 0.5*(ualpha(k-1) + ualpha(k)), 0.5*(ubeta(k-1)+ubeta(k)), muS, data);
%         K3 = EoMpolarAD(time(k)+ hhh/2, xx + hhh/2*K2, 0.5*(T(k-1) + T(k)), 0.5*(ualpha(k-1) + ualpha(k)), 0.5*(ubeta(k-1)+ubeta(k)),  muS, data);
%         K4 = EoMpolarAD(time(k)+ hhh, xx + hhh*K3, T(k-1), ualpha(k-1), ubeta(k-1), muS, data);    
% 
% %         Xprop(k+1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
%         Xprop(k-1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
%     end
%     
%     csi = zeros(1, (N-1)*7);
% 
%     for k = 1:N-1
%         csi((k-1)*7 + 1:(k-1)*7 + 7) = Xprop(k,:) - x(k,:);
%     end
    %collocation constraints
    DDD = zeros(1,(N-1)*7);
    
    %computation of collocation point
    for k = 1:N-1

        h = time(k+1) - time(k);

        xk   = x(k,:); 
        xkk  = x(k+1,:); 
        yk   = EoMpolarAD(time(k), xk, T(k), ualpha(k), ubeta(k), muS, data);
        ykk  = EoMpolarAD(time(k+1), xkk, T(k+1), ualpha(k+1), ubeta(k+1), muS, data);

        xc = 0.5*(xk + xkk) + h/8*(yk-ykk);    
        uc = ([T(k) ualpha(k) ubeta(k)] + [T(k+1) ualpha(k+1) ubeta(k+1)])/2;
        
        yc = EoMpolarAD(time(k+1), xc, uc(1), uc(2), uc(3), muS, data);

        DDD((k-1)*7 + 1:(k-1)*7 + 7) =  xk - xkk + h/6*(yk + 4*yc + ykk);
    
    end   
    
    % --- T modulation over distance, time (degradation), angle wrt sun
    switch data.ThrustModulationFlag
        case 0
            eta = 1;
        case 1
            % distance
            eta_space = solar_irradiance_cooler(x(:,1), 1, 1); % r is in AU
            eta = eta_space/max(eta_space);
        case 2
            % distance
            eta_space = solar_irradiance_cooler(x(:,1), 1, 1); % r is in AU
            
            % degradation
            % REFERENCES:
            D = 0.01674;
            time_in_years = time*sim.TU/(86400*365);
            eta_degradation = (1-D).^time_in_years;

            %total efficency
            eta = eta_space.*eta_degradation;
            eta = eta/max(eta);
        case 3
            % distance
            eta_space = solar_irradiance_cooler(x(:,1), 1, 1); % r is in AU
            
            % degradation
            % REFERENCES:
            D = 0.01674;
            time_in_years = time*sim.TU/(86400*365);
            eta_degradation = (1-D).^time_in_years;

            % angle of panels respect to sun
            % --- local cartesian
            [xx,yy,zz] = pol2cart(x(:,2),x(:,1),x(:,3));
            in_plane_position = [xx,yy,zeros(N,1)];
            in_plane_position_norm = vecnorm(in_plane_position,2,2);
            in_plane_position_vers = in_plane_position./in_plane_position_norm;
        
            z_vers = [0,0,1].*ones(N,3);
            tg_vers = -cross(in_plane_position_vers,z_vers);
        
            % ----- cartesian 3D thrust
            % beta = elevation, alpha = thrust angle in plane
            % x -> T t , y = -T rad, z -> T outofplane
            [Tx, Ty, Tz] = sph2cart(ualpha',ubeta',T');
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
            TT = rotate_local2ecplitic(r_in,T_cart_local,N,Href);
        
            r_transf_orbit  = [x(:,1).*cos(x(:,2)), x(:,1).*sin(x(:,2)), x(:,3)];
            RR = rotate_local2ecplitic(r_in,r_transf_orbit,N,Href);
        
        	% RR,TT, must be Nx3, cartesian position and thrust
            a_sun = acos(dot(RR,TT,2)./(vecnorm(TT,2,2).*x(:,1))) -pi/2; 
            eta_angle = cos(a_sun);
            eta_angle(isnan(eta_angle)) = cosd(4);

            %total efficency
            eta = eta_space.*eta_angle.*eta_degradation;
            eta = eta/max(eta);
    end
    
    % -- Thrust modulation and checks
    T_max_modulated = data.Tmax*eta;
    
    T_res = T - T_max_modulated'; % T is row, because the residual vector is build with a really big row
    check_T_ok = T_res<0;
    for i=1:data.n_int
        if check_T_ok(i) == 1
            T_res(i) = 0;
        else
           T_res(i) = T_res(i);
        end
    end
    
    %quality constrataints
%     ceq     = [csi, DDD, BCi, BCf]; 
    ceq     = [DDD,T_res];%, BCi, BCf]; 

    c = [];
end
