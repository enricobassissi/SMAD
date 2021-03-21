clear
clc
addpath time
addpath function

%%
date_ed = [2020, 1, 1, 0, 0, 0];
date_ld =  [2025, 1, 1, 0, 0, 0];
date_ea =  [2035, 1, 1, 0, 0, 0];
date_la =  [2060, 1, 1, 0, 0, 0];

mjd2000_ed = date2mjd2000(date_ed);
mjd2000_ld = date2mjd2000(date_ld);
mjd2000_ea = date2mjd2000(date_ea);
mjd2000_la = date2mjd2000(date_la);

%%
A=[ 1, -1];
b= -365.25*8;
lb=[ mjd2000_ed, mjd2000_ea];
ub=[ mjd2000_ld, mjd2000_la];
fun=@ff_ME2NE_nopga;

rng default
options = optimoptions(@ga);
options.PopulationSize = 80;
options.MaxGenerations = 150;
options.MaxStallGenerations = 60;
options.FunctionTolerance = 1e-6;
options.ConstraintTolerance = 1e-7;
options.MaxStallTime = 40;
options.MaxTime = 90;
options.UseParallel = true;
options.MigrationDirection = 'both';
options.MutationFcn = @mutationadaptfeasible;
options.CreationFcn = @gacreationnonlinearfeasible;
options.PlotFcn = {@gaplotbestf,@gaplotstopping};

fminconOptions = optimoptions(@fmincon,'PlotFcn',{'optimplotfval','optimplotx'},...
    'Display','iter','Algorithm', 'sqp','TolFun',1e-10);
options = optimoptions(options,'HybridFcn',{@fmincon, fminconOptions});

[T,dv,exitflag,output,population,scores]=ga(fun,2,A,b,[],[],lb,ub,[],options);

%%
min_dv=dv;
T_min=T;
date_start = mjd20002date(T_min(1));
date_arrival= mjd20002date(T_min(2));

fprintf ('\n start date from mercury is [%g %g %g %g %g %g] .\n',...
    date_start)
fprintf ('\n arrival date on neptune is [%g %g %g %g %g %g] .\n',...
    date_arrival)

fprintf ('\n the minimum dv is [%g] km/s.\n',...
    min_dv)

t1_sec=T_min(1)*24*3600;
t2_sec=T_min(2)*24*3600;

%% plot arcs

[kep_ME,ksun] = uplanet(T_min(1), 1);
[rm_ME, v1] = sv_from_coe(kep_ME,ksun);


[kep_NE,ksun] = uplanet(T_min(2), 8);
[rm_NE, v2] = sv_from_coe(kep_NE,ksun);

ToF12_lamb_min=t2_sec-t1_sec;
[a12_min,P,E,ERROR,VI12_min,VF12_min,TPAR,THETA] = lambertMR(rm_ME,rm_NE,ToF12_lamb_min,ksun,0,0,0,0);

% FULL ORBIT IN THE YEAR OF THE ACTUAL ENCOUNTER
mjd2000_start_ME = date2mjd2000(date_start);
% one year of Neptune in days
oneyearNE=165*365;
%PLOT MARS FULL ORBIT 
mjd2000_on_MA = date2mjd2000([2025, 1, 1, 0, 0, 0]);
% one year of Mars
oneyearMA=687;

%PLOT MERCURY FULL ORBIT 
mjd2000_on_NE = date2mjd2000(date_arrival );
% one year of Mercury
oneyearME=88;
% for full orbit propagation
n=100;
T_NEP=linspace(mjd2000_on_NE,mjd2000_on_NE+oneyearNE,n);

T_MARS=linspace(mjd2000_on_MA,mjd2000_on_MA+oneyearMA,n);

T_MERC=linspace(mjd2000_start_ME,mjd2000_start_ME+oneyearME,n);

% INITIAL AND FINAL ORBIT DEFINITION
% ID NEPTUNE 8, ID MARS 4, ID MERCURY 1
% [kep,ksun] = uplanet(mjd2000, ibody)

% COMPLETE ORBITS
for k=1:n
    [kep_NEP,ksun] = uplanet(T_NEP(k), 8);
    [r_NEP, v_NEP] = sv_from_coe(kep_NEP,ksun);  
    R_NEP(k,:)=r_NEP;
    V_NEP(k,:)=v_NEP;
    [kep_MARS,ksun] = uplanet(T_MARS(k), 4);
    [r_MARS, v_MARS] = sv_from_coe(kep_MARS,ksun);  
    R_MARS(k,:)=r_MARS;
    V_MARS(k,:)=v_MARS;
    [kep_MERC,ksun] = uplanet(T_MERC(k), 1);
    [r_MERC, v_MERC] = sv_from_coe(kep_MERC,ksun);  
    R_MERC(k,:)=r_MERC;
    V_MERC(k,:)=v_MERC;
end

% PLOT ORBITS AND BEST LAMBERT TRANSFER 
figure()

AU = astroConstants(2);

hold on
plot3(R_NEP(:,1)./AU,R_NEP(:,2)./AU,R_NEP(:,3)./AU,'--b');
plot3(R_MARS(:,1)./AU,R_MARS(:,2)./AU,R_MARS(:,3)./AU,'--r');
plot3(R_MERC(:,1)./AU,R_MERC(:,2)./AU,R_MERC(:,3)./AU,'--k');

t012 = t1_sec;
tf12 = t2_sec;
y012 = [rm_ME, VI12_min]';
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t12,y12] = ode113(@rates, [t012 tf12], y012, options, 'sun');
plot3( y12(:,1)./AU, y12(:,2)./AU, y12(:,3)./AU,'Color','g')

plot3(0,0,0,'*','Color','y')
plot3(rm_ME(1)./AU,rm_ME(2)./AU,rm_ME(3)./AU,'o','Color','k','MarkerSize',4);
plot3(rm_NE(1)./AU,rm_NE(2)./AU,rm_NE(3)./AU,'o','Color','b');
axis equal
grid minor
%legend ('full neptune','full mars','full mercury','trasnfer ME-MA','transfer MA-NE','sun','Location','southeast')
xlabel('AU')
ylabel('AU')
zlabel('AU')