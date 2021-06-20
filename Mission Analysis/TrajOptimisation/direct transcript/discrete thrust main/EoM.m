function [c, ceq] = EoM(X, data,sim)
    muS = data.muS;
    time = data.time;
    N = data.n_int;
        
    x = zeros(N,7);
    m = zeros(1, N);
    T = m; ualpha = m; ubeta = m;
    EoMpolarAD=@(t_hand, x_hand, T_hand, alpha_hand, beta_hand, muS_hand, data_hand) EoMpolar(t_hand, x_hand, T_hand, alpha_hand, beta_hand, muS_hand, data_hand,sim);
    for k=1:N
        %adimensional variables
        x(k,:)      = X((k-1)*10 + 1:(k-1)*10+7);
        % dimensional variables
        if X((k-1)*10+8)>=0.5
            T(k)=0.020;
        else
            T(k)=0;
        end
%         T(k)        = X((k-1)*10+8); 
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
    
    %quality constrataints
%     ceq     = [csi, DDD, BCi, BCf]; 
    ceq     = [DDD];%, BCi, BCf]; 

    c = [];
end
