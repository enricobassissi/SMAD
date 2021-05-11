%% --------------------------------------------------------------------- %%
%% --------------------------------- MARGINS --------------------------- %%
%% --------------------------------------------------------------------- %%
%% Paths
%% add path of functions and python stuff
str_path=split(pwd, 'TrajOptimisation\Impulsive\main');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'main');
imp_path=string(str_path(1));
addpath(genpath(imp_path));


%% Description
%{
Margins specified in this document are to be calculated using the following definition (unless
specifically defined differently):
margin% = (V2-V1)/V1
where V1 represents the initial value of the parameter without margin, 
and V2 the resulting value with margin. 
%}

% Use case of flyby 1+4 config
% load previous optimisation decided
load('soo_ps_flyby_3-33_1dayeachpointTraj.mat')

%% MAR-DV-010 
%{
5% or 10 m/s delta-v margin (whichever is highest) shall be applied to
deterministic, accurately calculated manoeuvres (trajectory manoeuvres as well
as detailed orbit maintenance manoeuvres) documented in a MAG/CReMA.
Note: 10 m/s are to be compared to 5% of the sum of all deterministic delta-v
terms and added only once if higher. 
%}
dV_not_margined = [sol.dV_extra_launch,sol.dVast1,sol.dVast2,sol.dVast3];
check_if_bigger_than_10ms = dV_not_margined*0.05;
% dV_margin_DV_010 = 0;
% 
% for i=1:length(dV_not_margined)
%     if dV_not_margined(i)*0.05 < 0.010 % 10 m/s
%         dV_margin_DV_010 = dV_margin_DV_010 + 0.01;
%     else
%         dV_margin_DV_010 = dV_margin_DV_010 + dV_not_margined(i)*0.05;
%     end
% end
dV_margin_DV_010 = sum(check_if_bigger_than_10ms);

%% MAR-DV-060 
%{ 
Gravity losses shall be quantified and added to the specified effective delta-v
when applicable (for instance, impulsive manoeuvres performed by chemical
propulsion engines).
Note 1: For gravity losses of deterministic nature, i.e. calculated based on the
thrust-to-mass ratio, MAR-DV-010 shall be applied.
Note 2: For estimated gravity losses where no thrust-to-mass ratio is available, the
margin on the gravity losses can be included in the estimate but shall be explicitly
stated. 
%}

%% MAR-DV-080 
%{
Based on the type of mission, a launcher dispersion correction manoeuvre
shall be included in the delta-v budget as follows:
- Any Earth orbit (from LEO to HEO): An allocation shall be estimated
based on the selected launchers’ injection accuracy (at 2 or 3 sigma
precision as agreed with the study team and over the complete launch
window)
- Direct escape launch: An allocation of 30 m/s shall be considered unless
more precise information is available.
%}
dV_margin_DV_080 = 0.030; % km/s

%% MAR-DV-090 
%{
In case of gravity assist manoeuvres (GAM), which are of stochastic nature, an
allocation of either 15 m/s (planetary GAMs) or 10 m/s (GAMs of planetary
moons) shall be added to the delta-v budget for each GAM performed by chemical
propulsion in order to account for preparation and correction of these
manoeuvres.
Note: These numbers are based on flight performance of Rosetta, previous SolO
analysis, and assessment made for JUICE for Venus Gravity Assist Manoeuvres.
Should the applicability be in doubt, these numbers shall be verified for the
relevant mission profile. 
%}
dV_margin_DV_090 = 0.015; % km/s

%% MAR-DV-100 
%{
In case of interplanetary approach navigation manoeuvres (towards the
Earth, a planet or a Moon) which are of stochastic nature, an allocation of 10 m/s
shall be added to the delta-v budget, if chemical propulsion is used.
Note: The interplanetary approach manoeuvre refers to the manoeuvre required
for the planetary approach that leads to arrival or orbit insertion.
%}
dV_margin_DV_100 = 0.010; % km/s

%% MARGINED
dV_margined = sum([sol.dV_extra_launch,sol.dVast1,sol.dVast2,sol.dVast3])+dV_margin_DV_010+dV_margin_DV_080

%% MAR-CP-010 
%{
The volume of the propellant tanks shall be sized for the total propellant
mass (including margins as defined in MAR-MAS-080) plus at least 10%.
Note: The result of this shall stay below the maximum qualified filling ratio. 
%}

%% MAR-MAS-010 
%{
Launcher margins (if any) are mission dependant and shall be defined on a
case by case basis and agreed with the study team. Launcher margins (if any)
shall be explicitly visible and shall be applied on the actual launcher vehicle
performance after subtraction of the launch adapter mass. Guidelines are:
- If the performance in the exact injection orbit is directly quoted in the
user manual, the launcher margin is not required.
- If the performance is calculated for a specific orbit or trajectory (e.g. for
interplanetary missions, as a function of the C3, declination and launch
date) by anyone else than the launcher authority (e.g. ESOC), a launcher
margin of 2% (5% maximum) should be included. 
%}

%% MAR-MAS-080 
%{
The propellant mass shall include 2% for propellant residuals.
Note 1: Propellant residuals cannot be used and account for left-over propellant
at EoL (e.g. propellant remaining in the fuel lines and at the bottom of the
tanks). Hence, residuals are accelerated during each manoeuvre during the
complete mission lifetime of the corresponding spacecraft.
Note 2: This margin is to be accounted by the expert assessing the propellant
mass.
Note 3: Propellant residuals are part of the total propellant mass which
determines the propellant tank sizing. 
%}

%% FINAL
dV_margined = sum(dV_not_margined)+dV_margin_DV_010+dV_margin_DV_080;