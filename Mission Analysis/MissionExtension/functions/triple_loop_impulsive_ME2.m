function [dvtot,time_vec_contour,time_to_go_vect,base_time_vec,horizons_best_ast_data,...
    TL_last_ast] = triple_loop_impulsive_ME2(sim,N,...
    end_of_mission,last_ast_name,names_ast_best_performing)

% --- TRIPLE LOOP METHOD, OVER ASTEROID & DEPARTURE & TOF --- %

% --- POSITIONS
% Asteroids
mjd2000_start = pystr2mjd2000(end_of_mission);
% for i=1:length(date_of_ca)
%     temp_max_date_ca(i) = pystr2mjd2000(date_of_ca(i));
% end
% mjd2000_stop = 1.2*max(temp_max_date_ca);
mjd2000_stop = mjd2000_start + 3*N;
epoch_stop = mjd20002pystr(mjd2000_stop);
PointOfView = 'Sun';
step = '2d';
type_elements = 'Vectors';
base_time_vec = [mjd2000_start:2:mjd2000_stop]';

DIS = length(base_time_vec); % discretization

% last_ast
py_TL_last_ast = py.neo_api_function.get_horizons_ephemerides(py.str(last_ast_name),py.str(PointOfView),...
                  py.str(end_of_mission),py.str(epoch_stop),py.str(step),py.str(type_elements));
TL_last_ast = double(py_TL_last_ast); % [x,y,z] in AU; [vx,vy,vz] in AU/day

% --- data extraction section of the other asteroids
% max_TOF = 1.2*(max(temp_max_date_ca) - mjd2000_start);
max_TOF = mjd2000_stop - mjd2000_start;
epoch_stop_other_ast_TL = mjd20002pystr(mjd2000_stop+max_TOF);

% complete orbit
for name=1:length(names_ast_best_performing)
    py_best_ast_data = py.neo_api_function.get_horizons_ephemerides(py.str(names_ast_best_performing(name)),py.str(PointOfView),...
                  py.str(end_of_mission),py.str(epoch_stop_other_ast_TL),py.str(step),py.str(type_elements));
    horizons_best_ast_data{name} = double(py_best_ast_data); % [x,y,z] in AU; [vx,vy,vz] in AU/day

%     figure()
%     plot3(horizons_best_ast_data{name}(:,1),horizons_best_ast_data{name}(:,2),horizons_best_ast_data{name}(:,3))
%     legend('names_ast_best_performing(name)')
end

DIS2 = length(horizons_best_ast_data{1});
time_to_go_vect = linspace(0, max_TOF, DIS2); % days -> triple loop stuff
% --- Actual triple loop
% departure_mjd2000_vec = linspace(mjd2000_start,mjd2000_end,DIS); % mjd2000
% tic
for k=1:length(names_ast_best_performing) % asteroid arrival selected
%     ast_to_go = names_ast_best_performing(k);
    TL_ast = horizons_best_ast_data{k};
    for i=1:DIS % departure variation
%         day_dep = base_time_vec(i);
%         states_ast_1 = interp1(base_time_vec,TL_last_ast,day_dep)'; % [x,y,z] in AU; [vx,vy,vz] in AU/day
        states_ast_1 = TL_last_ast(i,:); % [x,y,z] in AU; [vx,vy,vz] in AU/day
        states_ast_1(4:6) = states_ast_1(4:6)/86400*sim.TU; % it's adim now
        for j=1:DIS2 % tof variation
            tof = time_to_go_vect(j);
            states_ast_2 = TL_ast(j,:);
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
for i=1:DIS
    time_vec_contour(i) = datenum(datetime(mjd20002date(base_time_vec(i))));
end
clearvars i

end