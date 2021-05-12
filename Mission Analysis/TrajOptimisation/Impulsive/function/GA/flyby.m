function delta_v_p = flyby(Radius_of_planet, mu_flyby_planet, R_lim_from_planet, flyby_mjd2000, v_minus, v_plus)
[kep,ksun] = uplanet(flyby_mjd2000,3);
            [~,V]=sv_from_coe(kep,ksun);
            v_inf_before=v_minus-V';
            v_inf_after=v_plus-V';
            delta=acos(dot(v_inf_before,v_inf_after)/(norm(v_inf_before)*norm(v_inf_after)));
%             eq =  @(x) asin(1/(1 + (x^2+flyby_planet.R+440)*norm(v_inf_before)^2/flyby_planet.mu)) + asin(1/(1 + (x^2+flyby_planet.R+440)*norm(v_inf_after)^2/flyby_planet.mu)) - delta;
%             rp_ = fzero(eq,1);
%             rp=(rp_^2+flyby_planet.R+440);
%             if isnan(rp)
%                 delta_v_p=10^6;
%             else
            opt=optimset('Display', 'off');
%             opt=optimoptions('fsolve', 'TolFun', 1e-13, 'TolX', 1e-13,'Display','off');
            eq =  @(x) asin(1/(1 + x*norm(v_inf_before)^2/mu_flyby_planet)) + asin(1/(1 + x*norm(v_inf_after)^2/mu_flyby_planet)) - delta;
%             try
%                 rp = vpa_solve(eq,Radius_of_planet+441,opt);
                rp = fzero(eq,Radius_of_planet+R_lim_from_planet+1,opt);
%             catch
%                 rp = NaN;
%             end
            
            if or((rp<(Radius_of_planet+R_lim_from_planet)),isnan(rp))
                delta_v_p='Not found';
            else
            e_before=1+rp*norm(v_inf_before)^2/mu_flyby_planet;
%             delta_before=2*asin(1/e_before);
            a_before=rp/(1-e_before);
            e_after=1+rp*norm(v_inf_after)^2/mu_flyby_planet;
%             delta_after=2*asin(1/e_after);
            a_after=rp/(1-e_after);
            

%             switch 'try'
%                 case 'try'
%                     u=-cross(v_inf_before,v_inf_after)./norm(cross(v_inf_before,v_inf_after));
%                 case 'front'
%                     u=[0;0;1];
% %                 case 'behind'
%                     u=[0;0;-1];
% %                 case 'under'
%                     u=cross( [0;0;1], v_inf_before/norm(v_inf_before));
% %                 case 'over'
%                     u=cross( [0;0;-1], v_inf_before/norm(v_inf_before));
%             end
%             rot_vector=delta_before*(u'/norm(u'));
%             rot_matrix=rotationVectorToMatrix(rot_vector);
%             v_inf_after_noprop=rot_matrix*v_inf_before;
%             delta_v_noprop=v_inf_after_noprop-v_inf_before;
%             delta_v_versor_noprop=delta_v_noprop./(norm(delta_v_noprop));

            vp_norm_before=sqrt(mu_flyby_planet*((2/rp)-(1/a_before)));
            vp_norm_after=sqrt(mu_flyby_planet*((2/rp)-(1/a_after)));
            delta_v_p=abs(vp_norm_after-vp_norm_before);
            %for plotting etc
%             if(bool_plot)
%             vp_before=cross(-delta_v_versor_noprop,u).*(vp_norm_before);
%             vp_after=cross(-delta_v_versor_noprop,u).*(vp_norm_after);
%             given_delta_v=v_inf_after-v_inf_before;
%             T_final=5000;
%             figure(2)
%             propagation=obj.propagate_around(flyby_planet, -rp*delta_v_versor_noprop, vp_after, T_final);
%             propagation_bw=obj.propagate_backward_around(flyby_planet, -rp*delta_v_versor_noprop, vp_before, T_final);
%             propagation.plot_hp();
%             hold on
%             propagation_bw.plot_hp();
            end
end