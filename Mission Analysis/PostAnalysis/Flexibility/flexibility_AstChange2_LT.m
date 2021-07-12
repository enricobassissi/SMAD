%% --------------------------------------------------------------------- %%
%% ---------------------- FLEXIBILITY ANALYSIS ------------------------- %%
%% --------------------- Let's change one asteroid --------------------- %%
%% --------------------------------------------------------------------- %%
%% Description
% We optimise the solution for a given set of asteroids. But then after
% some observations of them we understand that one is not suitable for the
% mission at all. so we want to change it with one of the remaining
% asteroids in the list. How this changing affect the mission? We change
% one by one all the asteroids (keeping the sequence unique of course) and
% see the effect on overall dV and TOF.

%% Setup for default options
set(0, 'DefaultTextFontSize', 20)
set(0, 'DefaultAxesFontSize', 20)
set(0, 'DefaultLegendFontSize', 20)
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLineLineWidth', 1.8)
format short

%% Initializing the Environment
clear; close all; clc;

% Palette ESA
colors = [0    50   71;... % (1) DEEP SPACE
          207  29   57;... % (2) EXCITE RED +1
          0    142  122;... % (3) PURE TEAL +1
          251  171  24;... % (4) ENLIGHT YELLOW
          244  121  32;... % (5) ENLIGHT YELLOW +1
          150  1    54;... % (6) EXCITE RED +2
          167  85   52;... % (7) ENLIGHT YELLOW +2
          0    97   158;... % (8) TRUSTY AZURE +1
          30   51   120;... % (9) TRUSTY AZURE +2
          0    103  98;... % (10) PURE TEAL +2
          51   94   111;... % (11) DEEP SPACE -1
          0    174  157;... % (12) PURE TEAL
          0    0    0]./255; % (13) BLACK

%% add path of functions and python stuff
str_path=split(pwd, 'PostAnalysis\Flexibility');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'Flexibility');
imp_path=string(str_path(1));
addpath(genpath(imp_path));

%% Call to NASA JPL Horizons to get Asteroid's Ephemerides
% Import module of Python
try 
    module = py.importlib.import_module('neo_api_function');
catch
    copyfile(py_path+'\neo_api_function.py', pwd, 'f'); 
    module = py.importlib.import_module('neo_api_function');
end

%% Asteroids
load('data_elements_matrix_44_63_2SC.mat')
load('ws_2RL_all_indietro_moo2.mat')

%% simulation parameters
sim.mu_dim    = 132712440018              ; % actractor parameter [km^3 s^-2]
sim.DU        = 149597870.7               ; % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu        = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol     = 100;                    % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 

sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = -1;                     % direction of integration (1 FW, -1 BW), 
                                       % 1 is like imposing wet mass at beginning
sim.TOF_imposed_flag = 1;
sim.PS.Isp = 3200/sim.TU;  % non-dimensional specific impulse
% sim.PS.Isp = 4500/sim.TU;  % non-dimensional specific impulse % simone
sim.M1_end = 160; % SC wet mass [kg] %%
sim.M2_end = 160; % SC wet mass [kg] %%
sim.M_pods = 3.5; % mass of the payloads + landing stuff [kg] %%
sim.max_Available_Thrust = 0.02; % 5 [mN], BepiColombo is 250 mN but it's much bigger

%% Flexibility Parameters
ast_solution = [sol.asteroid_1,sol.asteroid_2,sol.asteroid_a,sol.asteroid_b];
ast_to_avoid = [];
ast_to_keep = [];
% Building the permutation matrix of the asteroid sequence but changing one
% at a time the 1st,2nd,ath,bth
for i=1:2 % rows, ID number of the SC
    for j=1:2 % columns, number of ast per each SC
        if i == 1 && j == 1 % first asteroid of the two
            ast_to_keep = ast_solution(1);
            ast_to_avoid = [ast_solution(3:4), ast_solution(2)];
            data.OtherSCAstFixedFlex{i,j} = ast_solution(3:4);
        elseif i == 1 && j == 2 % second asteroid of the two
            ast_to_keep = ast_solution(2);
            ast_to_avoid = [ast_solution(3:4), ast_solution(1)];
            data.OtherSCAstFixedFlex{i,j} = ast_solution(3:4);
        elseif i == 2 && j == 1 % first asteroid of the two
            ast_to_keep = ast_solution(3);
            ast_to_avoid = [ast_solution(1:2), ast_solution(4)];
            data.OtherSCAstFixedFlex{i,j} = ast_solution(1:2);
        elseif i == 2 && j == 2 % first asteroid of the two
            ast_to_keep = ast_solution(4);
            ast_to_avoid = [ast_solution(1:2), ast_solution(3)];
            data.OtherSCAstFixedFlex{i,j} = ast_solution(1:2);
        end
        data.AstToKeep{i,j} = ast_to_keep;
        TF_same = ast_to_avoid(3) == data.asteroid_names;
        TF1 = ast_to_avoid(1) == data.asteroid_names;
        TF2 = ast_to_avoid(2) == data.asteroid_names;
        TF = TF_same+TF1+TF2;
        data_elements_matrix_2SC = data.data_elements_matrix(~TF,:);

         all_perm_rep = [repelem(ast_to_keep,length(data_elements_matrix_2SC),1),data_elements_matrix_2SC(:,1);
                                           data_elements_matrix_2SC(:,1),repelem(ast_to_keep,length(data_elements_matrix_2SC),1)];
        [~,data.AstFlexSequenceMatrix{i,j},~] = sequences_local_pruning3(data_elements_matrix_2SC,data.p_number,all_perm_rep);     
    end
end
clearvars i j

%% Boundaries
% Departure dates (1)
bound_flex.date_ed = [2025, 1, 1, 0, 0, 0];
bound_flex.date_ld =  [2028, 1, 1, 0, 0, 0];
bound_flex.mjd2000_ed = date2mjd2000(bound_flex.date_ed)*3600*24/sim.TU;
bound_flex.mjd2000_ld = date2mjd2000(bound_flex.date_ld)*3600*24/sim.TU;
% TOF1 (2)
bound_flex.TOF1_min = 1*365*3600*24/sim.TU; %1*365
bound_flex.TOF1_max = 4*365*3600*24/sim.TU; %3*365
% TOF2 (3)
bound_flex.TOF2_min = 1*365*3600*24/sim.TU; %1*365
bound_flex.TOF2_max = 4*365*3600*24/sim.TU; %3*365
% TOFa (4)
bound_flex.TOFa_min = 1*365*3600*24/sim.TU; %1*365
bound_flex.TOFa_max = 4*365*3600*24/sim.TU; %3*365
% TOFb (5)
bound_flex.TOFb_min = 1*365*3600*24/sim.TU; %1*365
bound_flex.TOFb_max = 4*365*3600*24/sim.TU; %3*365
% N REV 1 (6)
bound_flex.N_REV1_min = 1; %0
bound_flex.N_REV1_max = 2; %3
% N REV 2 (7)
bound_flex.N_REV2_min = 1; %0
bound_flex.N_REV2_max = 2; %3
% N REV a (8)
bound_flex.N_REVa_min = 1; %0
bound_flex.N_REVa_max = 2; %3
% N REV b (9)
bound_flex.N_REVb_min = 1; %0
bound_flex.N_REVb_max = 2; %3
% C3 stuff
% Constraint on C3 Launcher (10)
sim.C3_max = 30; % km^2/s^2
bound_flex.v_inf_magn_min = 0;
bound_flex.v_inf_magn_max = sqrt(sim.C3_max)/sim.DU*sim.TU;
% azimuth (11)
bound_flex.az_min = -pi;
bound_flex.az_max = pi;
% elevation (12)
bound_flex.el_min = -pi/5;
bound_flex.el_max = pi/5;

% Coasting time 1 (13)
bound_flex.CT1_min = 30*3600*24/sim.TU; % 30 days
bound_flex.CT1_max = 80*3600*24/sim.TU; 
% Coasting time 2 (14)
bound_flex.CTa_min = 30*3600*24/sim.TU; % 30 days
bound_flex.CTa_max = 80*3600*24/sim.TU; 

% x = [MJD0,TOF1,TOF2,TOFa,TOFb,NREV,NREV2,NREVa,NREVb,IDP,v_inf_magn,az,el,CT1,CTa]
bound_flex.lb = [bound_flex.mjd2000_ed, bound_flex.TOF1_min, bound_flex.TOF2_min, ...
    bound_flex.TOFa_min, bound_flex.TOFb_min, bound_flex.N_REV1_min, bound_flex.N_REV2_min, ...
    bound_flex.N_REVa_min, bound_flex.N_REVb_min, ...
    bound_flex.v_inf_magn_min, bound_flex.az_min, bound_flex.el_min, bound_flex.CT1_min, bound_flex.CTa_min]; % Lower bound
bound_flex.ub = [bound_flex.mjd2000_ld, bound_flex.TOF1_max, bound_flex.TOF2_max, ...
    bound_flex.TOFa_max, bound_flex.TOFb_max, bound_flex.N_REV1_max, bound_flex.N_REV2_max, ...
    bound_flex.N_REVa_max, bound_flex.N_REVb_max, ...
    bound_flex.v_inf_magn_max, bound_flex.az_max, bound_flex.el_max, bound_flex.CT1_max, bound_flex.CTa_max]; % Upper bound

%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints

%% Options
options = optimoptions(@gamultiobj);
% options.PlotFcn = @gaplotpareto;
options.Display = 'none';
% y = score; -> phenotype, Measure the distance in fitness function space; 
% y = pop; -> genotype, Measure the distance in decision variable space.
options.DistanceMeasureFcn = {@distancecrowding,'phenotype'};
% A hybrid function is another minimization function that runs after the 
% multiobjective genetic algorithm terminates
% options.HybridFcn = 'fgoalattain';

% ---- remember to change the integer position in the vector
options.CreationFcn = @int_pop_2RL_moo_AC2_flex;
options.MutationFcn = @int_mutation_2RL_moo_AC2_flex;
options.CrossoverFcn = @int_crossoverarithmetic_2RL_moo_AC2_flex;

options.PopulationSize = 300; % ideal 1000
options.ParetoFraction = 0.6;
options.MaxGenerations = 50; % ideal 100

options.FunctionTolerance = 1e-9;
options.MaxStallGenerations = ceil(options.MaxGenerations/2);

% Parallel pool
% Open the parallel pool
par_pool = gcp; 
if isempty(par_pool)
    poolsize = 0;
else
    poolsize = par_pool.NumWorkers;
end

options.UseParallel = true;

numberOfVariables = length(bound_flex.ub); % Number of decision variables

%% Multi Objective Loop
% initialize
idx_i = 0; idx_j = 0; idx_k = 0;

for i=1:2 % rows, which spacecraft we are considering
    idx_i = idx_i+1;
    idx_j = 0;
    for j=1:2 % columns, which asteroids that that sc is visiting
        idx_j = idx_j+1;
        idx_k = 0;
        data.IDX_WhichSpacecraftWeAreChanging = i;
        data.IDX_WhichAstWeAreChanging = j;
        data.OtherSpacecraftFixedAst = data.OtherSCAstFixedFlex{i,j};
        data.PermutCycle = data.AstFlexSequenceMatrix{i,j};
        TFF = unique(data.PermutCycle);
        AstListToCheck = TFF(~(unique(data.PermutCycle)==data.AstToKeep{i,j}));
        data.AstListToCheck{i,j} = AstListToCheck;
        for k=1:length(AstListToCheck)
            idx_k = idx_k+1;
            fprintf('Cycle i: %d, Cycle j: %d, Cycle k: %d  \n',idx_i, idx_j, idx_k);
            
            data.NextAstSameSpacecraft = AstListToCheck(k);
            FitnessFunction = @(x) ff_2RL_all_indietro_moo_AC2_flex(x, data, sim); % Function handle to the fitness function
            [xx_flex,Fval_flex,exitFlag,Output,population,score] = gamultiobj(FitnessFunction,numberOfVariables,constr.A, ...
                constr.b,constr.Aeq,constr.beq,bound_flex.lb,bound_flex.ub,constr.nonlcon,options);
            FVAL_stored{i,j,k} = Fval_flex;
            XX_stored{i,j,k} = xx_flex;

            % find knee solution
            knee_sol_Fval = sqrt(Fval_flex(:,1).^2+(Fval_flex(:,2)./500).^2);
            idx_knee = find(min(knee_sol_Fval) == knee_sol_Fval);
            idx_knee = idx_knee(1);
            x_flex = xx_flex(idx_knee,:);
            X_stored{i,j,k} = x_flex;

            % save the result
            if data.IDX_WhichSpacecraftWeAreChanging == 1
                if data.IDX_WhichAstWeAreChanging == 1
                    % 1st sc ast change
                    asteroid_1 = data.AstToKeep{1,1};
                    asteroid_2 = data.NextAstSameSpacecraft;
                elseif data.IDX_WhichAstWeAreChanging == 2
                    % 1st sc ast change
                    asteroid_1 = data.NextAstSameSpacecraft;
                    asteroid_2 = data.AstToKeep{1,2};
                end
                % 2nd sc's asteroid are defined
                asteroid_a = data.OtherSpacecraftFixedAst(1);
                asteroid_b = data.OtherSpacecraftFixedAst(2);

            elseif data.IDX_WhichSpacecraftWeAreChanging == 2
                if data.IDX_WhichAstWeAreChanging == 1
                    % 2nd sc's asteroid are changing
                    asteroid_a = data.AstToKeep{2,1};
                    asteroid_b = data.NextAstSameSpacecraft;
                elseif data.IDX_WhichAstWeAreChanging == 2
                    % 2nd sc's asteroid are changing
                    asteroid_a = data.NextAstSameSpacecraft;
                    asteroid_b = data.AstToKeep{2,2};
                end
                % 1st sc ast fixed
                asteroid_1 = data.OtherSpacecraftFixedAst(1);
                asteroid_2 = data.OtherSpacecraftFixedAst(2);
            end

            flex_order{i,j,k}.asteroid_sequence_SC1 = [asteroid_1,asteroid_2];
            flex_order{i,j,k}.asteroid_sequence_SC2 = [asteroid_a,asteroid_b];
            flex_order{i,j,k}.ast_1 = flex_order{i,j,k}.asteroid_sequence_SC1(1);
            flex_order{i,j,k}.ast_2 = flex_order{i,j,k}.asteroid_sequence_SC1(2);
            flex_order{i,j,k}.ast_a = flex_order{i,j,k}.asteroid_sequence_SC2(1);
            flex_order{i,j,k}.ast_b = flex_order{i,j,k}.asteroid_sequence_SC2(2);
            flex_order{i,j,k}.MJD0 = x_flex(1);
            flex_order{i,j,k}.TOF1 = x_flex(2);
            flex_order{i,j,k}.TOF2 = x_flex(3);
            flex_order{i,j,k}.TOFa = x_flex(4);
            flex_order{i,j,k}.TOFb = x_flex(5);
            flex_order{i,j,k}.CT1 = x_flex(13);
            flex_order{i,j,k}.CTa = x_flex(14);
            flex_order{i,j,k}.dep_date = mjd20002date(x_flex(1))';
            flex_order{i,j,k}.TOF_tot_D_SC1 = x_flex(2)+x_flex(3)+x_flex(13);
            flex_order{i,j,k}.TOF_tot_Y_SC1 = flex_order{i,j,k}.TOF_tot_D_SC1/365;
            flex_order{i,j,k}.TOF_tot_D_SC2 = x_flex(4)+x_flex(5)+x_flex(14);
            flex_order{i,j,k}.TOF_tot_Y_SC2 = flex_order{i,j,k}.TOF_tot_D_SC2/365;
            flex_order{i,j,k}.end_of_mission_date_SC1 = mjd20002date(x_flex(1)+flex_order{i,j,k}.TOF_tot_D_SC1)';
            flex_order{i,j,k}.end_of_mission_date_SC2 = mjd20002date(x_flex(1)+flex_order{i,j,k}.TOF_tot_D_SC2)';

            flex_order{i,j,k}.Nrev = [x_flex(6),x_flex(7),x_flex(8),x_flex(9)];
            flex_order{i,j,k}.v_inf_magn = x_flex(11)*sim.DU/sim.TU;
            flex_order{i,j,k}.az = x_flex(11);
            flex_order{i,j,k}.el = x_flex(12);
            flex_order{i,j,k}.az_deg = rad2deg(x_flex(11));
            flex_order{i,j,k}.el_deg = rad2deg(x_flex(12));
            flex_order{i,j,k}.obj_fun1 = Fval_flex(idx_knee,1);
            flex_order{i,j,k}.obj_fun2 = Fval_flex(idx_knee,2);

            plot_values = plot_ff_2RL_all_indietro_AC2_flex(x_flex,sim,data);
            flex_order{i,j,k}.max_T = max([max(plot_values.T_1_magn),max(plot_values.T_2_magn),max(plot_values.T_a_magn),max(plot_values.T_b_magn)]);
            flex_order{i,j,k}.max_MF = max(plot_values.mass_fract_SC2,plot_values.mass_fract_SC2);
        end
    end
end
el_time=toc;
clearvars i j k idx_i idx_j idx_k

%% plot the result
figure()
idx_color = 0;
for i=1:2 % rows
    for j=1:2 % columns
        KK = data.AstListToCheck{i,j};
        idx_color = idx_color+1;
        MF_flex_order = zeros(length(KK),1);
        T_flex_order = zeros(length(KK),1);
        
%         data.IDX_WhichSpacecraftWeAreChanging = i;
%         data.IDX_WhichAstWeAreChanging = j;
%         data.OtherSpacecraftFixedAst = data.OtherSCAstFixedFlex{i,j};
%         data.PermutCycle = data.AstFlexSequenceMatrix{i,j};
%         TFF = unique(data.PermutCycle);
%         AstListToCheck = TFF(~(unique(data.PermutCycle)==data.AstToKeep{i,j}));
%         data.AstListToCheck{i,j} = AstListToCheck;
        
        for k=1:length(KK)
%             data.NextAstSameSpacecraft = AstListToCheck(k);
%             [sol] = plot_ff_2RL_all_indietro_AC2_flex(X_stored{i,j,k},sim,data);
            MF_flex_order(k,1) = flex_order{i,j,k}.max_MF;
            T_flex_order(k,1) = flex_order{i,j,k}.max_T*1000;
%             MF_flex_order(k,1) = max(flex_order{i,j,k}.obj_fun1);
%             T_flex_order(k,1) = max(flex_order{i,j,k}.obj_fun2);
            if MF_flex_order(k,1) == 0 || MF_flex_order(k,1) > 0.6
                MF_flex_order(k,1) = NaN;
            end
            if T_flex_order(k,1) == 0 || T_flex_order(k,1) > 400
                T_flex_order(k,1) = NaN;
            end
        end
        scatter(MF_flex_order(:,1),T_flex_order(:,1),[],colors(idx_color,:),...
            'filled','DisplayName', num2str([i,j]));
        hold on; grid on;
        xlabel('MF [-]'); ylabel('Thrust [mN]'); 
        legend('show','Location','northoutside','NumColumns',2)

%         OBJ1(i,1) = flex{i}.max_MF;
%         if OBJ1(i,1) == 0 || OBJ1(i,1) > 0.4
%             OBJ1(i,1) = NaN;
%         end
    end
end
clearvars i j k KK idx_color


% figure('Name','Flexibility over Departure Time')
% yyaxis left
% plot(date_vect,OBJ1,'*','Color',colors(1,:))
% hold on
% plot(datenum(datetime(mjd20002date(sol.departure_mjd2000))),max(sol.mass_fract_SC2, sol.mass_fract_SC1),'*','Color',colors(3,:))
% ylabel('Mass Fraction')
% set(gca,'ycolor',colors(1,:)) 
% yyaxis right
% plot(date_vect,OBJ2,'o','Color',colors(2,:))
% plot(datenum(datetime(mjd20002date(sol.departure_mjd2000))),max(sol.max_T_SC1, sol.max_T_SC2)*1000,'o','Color',colors(3,:))
% ylabel('Max Thrust [mN]')
% set(gca,'ycolor',colors(2,:)) 
% % xlabel('Departure MJD2000')
% ax = gca;
% ax.XTick=date_vect(1:4:end) ;
% xtickangle(30)
% datetick('x','yyyy mmm dd','keepticks')
% xlim ([date_vect(1) date_vect(end)])
% xline (datenum(datetime(mjd20002date(sol.departure_mjd2000))),'--','Color',colors(3,:),'LineWidth',2)


% figure('Name','Changing and optimising the asteroid sequence order')
% scatter(MF_flex_order,T_flex_order,[],colors(1,:),'filled')
% scatter(MF_flex_order(:,1),T_flex_order(:,1),[],colors(1,:),'filled',...
%     'DisplayName','Ast 1 Variations and Reordering'); % all changes on first ast
% hold on
% scatter(MF_flex_order(:,2),T_flex_order(:,2),[],colors(2,:),'filled',...
%     'DisplayName','Ast 2 Variations and Reordering'); % all changes on 2nd ast
% scatter(MF_flex_order(:,3),T_flex_order(:,3),[],colors(4,:),'filled',...
%     'DisplayName','Ast 3 Variations and Reordering'); % all changes on 3rd ast
% scatter(MF_flex_order(:,4),T_flex_order(:,4),[],colors(12,:),'filled',...
%     'DisplayName','Ast 4 Variations and Reordering'); % all changes on 4th ast
% xlabel('$MF$ [-]'); ylabel('Thrust [N]'); 
% legend('show','Location','northoutside','NumColumns',2)