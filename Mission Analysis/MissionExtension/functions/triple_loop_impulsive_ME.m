function [dvtot,time_vec_contour,time_to_go_vect,base_time_vec,horizons_best_ast_data,...
    TL_last_ast] = triple_loop_impulsive_ME(sim,N,...
    end_of_mission,date_of_ca,last_ast_name,names_ast_best_performing)

% --- TRIPLE LOOP METHOD, OVER ASTEROID & DEPARTURE & TOF --- %

% --- POSITIONS
% Asteroids
mjd2000_start = pystr2mjd2000(end_of_mission);
for i=1:length(date_of_ca)
    temp_max_date_ca(i) = pystr2mjd2000(date_of_ca(i));
end
mjd2000_stop = 1.2*max(temp_max_date_ca);
epoch_stop = mjd20002pystr(mjd2000_stop);
PointOfView = 'Sun';
step = '2d';
type_elements = 'Vectors';
base_time_vec = [mjd2000_start:2:mjd2000_stop]';

% last_ast
py_TL_last_ast = py.neo_api_function.get_horizons_ephemerides(py.str(last_ast_name),py.str(PointOfView),...
                  py.str(end_of_mission),py.str(epoch_stop),py.str(step),py.str(type_elements));
TL_last_ast = double(py_TL_last_ast); % [x,y,z] in AU; [vx,vy,vz] in AU/day

% --- data extraction section
% complete orbit
for name=1:length(names_ast_best_performing)
    py_best_ast_data = py.neo_api_function.get_horizons_ephemerides(py.str(names_ast_best_performing(name)),py.str(PointOfView),...
                  py.str(end_of_mission),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_best_ast_data{name} = double(py_best_ast_data); % [x,y,z] in AU; [vx,vy,vz] in AU/day

%     figure()
%     plot3(horizons_best_ast_data{name}(:,1),horizons_best_ast_data{name}(:,2),horizons_best_ast_data{name}(:,3))
%     legend('names_ast_best_performing(name)')
end

% --- Actual triple loop
max_TOF = 1.2*(max(temp_max_date_ca) - mjd2000_start);
time_to_go_vect = linspace(0, max_TOF, 2*N); % days
departure_mjd2000_vec = linspace(mjd2000_start,1.01*max(temp_max_date_ca),N); % mjd2000
% tic
for k=1:length(names_ast_best_performing) % asteroid arrival selected
%     ast_to_go = names_ast_best_performing(k);
    TL_ast = horizons_best_ast_data{k};
    for i=1:length(departure_mjd2000_vec) % departure variation
        day_dep = departure_mjd2000_vec(i);
        states_ast_1 = interp1(base_time_vec,TL_last_ast,day_dep)'; % [x,y,z] in AU; [vx,vy,vz] in AU/day
        states_ast_1(4:6) = states_ast_1(4:6)/86400*sim.TU; % it's adim now
        for j=1:length(time_to_go_vect) % tof variation
            tof = time_to_go_vect(j);
            states_ast_2 = interp1(base_time_vec,TL_ast,day_dep+tof)';
            states_ast_2(4:6) = states_ast_2(4:6)/86400*sim.TU; % it's adim now
            [dv,~,~]=lambert_solver_rendezvous(states_ast_1(1:3),states_ast_2(1:3),...
                states_ast_1(4:6),states_ast_2(4:6),tof*86400/sim.TU,sim.mu); % everything adim... AU, adim vel, adim mu
            dvtot(i,j,k) = dv*sim.DU/sim.TU; % ri adimensionalise -> km/s
%             states_ast_2(4:6) = states_ast_2(4:6)/86400*sim.DU;
%             [dvtot(i,j,k),~,~]=lambert_solver_rendezvous(states_ast_1(1:3)*sim.DU,states_ast_2(1:3)*sim.DU,...
%                 states_ast_1(4:6),states_ast_2(4:6),tof*86400,sim.mu_dim);
        end
    end
end
% el_time = toc;
clearvars i j k
for i=1:length(departure_mjd2000_vec)
    time_vec_contour(i) = datenum(datetime(mjd20002date(departure_mjd2000_vec(i))));
end
clearvars i

end