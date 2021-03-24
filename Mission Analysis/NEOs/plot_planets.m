function [h, EA] = plot_planets(mjd2000_start_sim, varargin)

AU = astroConstants(2);
muSun = astroConstants(4);

oneyearME = 88; % one year of Mercury in days
oneyearVE = 225; % one year of Mercury in days
oneyearEA = 365; % one year of Earth in days
oneyearMA = 687; % one year of Mars in days
oneyearJU = 2*11.862*365; % one year of Jupiter in days
    
if isempty(varargin)
    % for full orbit propagation
    n = 100;
    T_ME = linspace(mjd2000_start_sim, mjd2000_start_sim + oneyearME, n);
    T_VE = linspace(mjd2000_start_sim, mjd2000_start_sim + oneyearVE, n);
    T_EA = linspace(mjd2000_start_sim, mjd2000_start_sim + oneyearEA, n);
    T_MA = linspace(mjd2000_start_sim, mjd2000_start_sim + oneyearMA, n);
    T_JU = linspace(mjd2000_start_sim, mjd2000_start_sim + oneyearJU, n*2);
elseif length(varargin) == 1
    mjd2000_end_sim = cell2mat(varargin(1));
    n = ceil((mjd2000_end_sim - mjd2000_start_sim)/10);
    T_ME = linspace(mjd2000_start_sim, mjd2000_end_sim, n);
    T_VE = linspace(mjd2000_start_sim, mjd2000_end_sim, n);
    T_EA = linspace(mjd2000_start_sim, mjd2000_end_sim, n);
    T_MA = linspace(mjd2000_start_sim, mjd2000_end_sim, n);
    T_JU = linspace(mjd2000_start_sim, mjd2000_end_sim, n*2);
elseif length(varargin) == 2
    mjd2000_end_sim = cell2mat(varargin(1));
    n = cell2mat(varargin(2));
    T_ME = linspace(mjd2000_start_sim, mjd2000_end_sim, n);
    T_VE = linspace(mjd2000_start_sim, mjd2000_end_sim, n);
    T_EA = linspace(mjd2000_start_sim, mjd2000_end_sim, n);
    T_MA = linspace(mjd2000_start_sim, mjd2000_end_sim, n);
    T_JU = linspace(mjd2000_start_sim, mjd2000_end_sim, n*2);
else
    disp('too many or no input argument')
end

R_ME = zeros(n, 3); V_ME = zeros(n, 3);
R_VE = zeros(n, 3); V_VE = zeros(n, 3);
R_EA = zeros(n, 3); V_EA = zeros(n, 3);
R_MA = zeros(n, 3); V_MA = zeros(n, 3);
R_JU = zeros(n, 3); V_JU = zeros(n, 3);

for k=1:n
    [kep_ME, muSun] = uplanet(T_ME(k), 1);
    [r_ME, v_ME] = sv_from_coe(kep_ME, muSun);  
    R_ME(k, :) = r_ME;
    V_ME(k, :) = v_ME;
    
    [kep_VE, muSun] = uplanet(T_VE(k), 2);
    [r_VE, v_VE] = sv_from_coe(kep_VE, muSun);  
    R_VE(k, :) = r_VE;
    V_VE(k, :) = v_VE;
    
    [kep_EA, muSun] = uplanet(T_EA(k), 3);
    [r_EA, v_EA] = sv_from_coe(kep_EA, muSun);  
    R_EA(k, :) = r_EA;
    V_EA(k, :) = v_EA;
    
    [kep_MA, muSun] = uplanet(T_MA(k), 4);
    [r_MA, v_MA] = sv_from_coe(kep_MA, muSun);  
    R_MA(k, :) = r_MA;
    V_MA(k, :) = v_MA;
    
    [kep_JU, muSun] = uplanet(T_JU(k), 5);
    [r_JU, v_JU] = sv_from_coe(kep_JU, muSun);  
    R_JU(k, :) = r_JU;
    V_JU(k, :) = v_JU;
end

% PLOT ORBITS
h.ME = plot3(R_ME(:,1)./AU, R_ME(:,2)./AU, R_ME(:,3)./AU, '--k');
hold on
h.VE = plot3(R_VE(:,1)./AU, R_VE(:,2)./AU, R_VE(:,3)./AU, '--m');
h.EA = plot3(R_EA(:,1)./AU, R_EA(:,2)./AU, R_EA(:,3)./AU, '--g');
h.MA = plot3(R_MA(:,1)./AU, R_MA(:,2)./AU, R_MA(:,3)./AU, '--r');
h.JU = plot3(R_JU(:,1)./AU, R_JU(:,2)./AU, R_JU(:,3)./AU, '--y');

EA.r = R_EA;

end