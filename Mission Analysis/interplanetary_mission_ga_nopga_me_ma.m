clear
clc
addpath time
addpath function

%%
date_ed = [2021, 1, 1, 0, 0, 0];
date_ld =  [2025, 1, 1, 0, 0, 0];
date_ea =  [2021, 1, 1, 0, 0, 0];
date_la =  [2030, 1, 1, 0, 0, 0];

mjd2000_ed = date2mjd2000(date_ed);
mjd2000_ld = date2mjd2000(date_ld);
mjd2000_ea = date2mjd2000(date_ea);
mjd2000_la = date2mjd2000(date_la);

%%
A=[ 1, -1];
b= -365.25*0.1;
lb=[ mjd2000_ed, mjd2000_ea];
ub=[ mjd2000_ld, mjd2000_la];
fun=@ff_ME2MA_nopga;

rng default
options = optimoptions(@ga);
options.PopulationSize = 80;
options.MaxGenerations = 150;
options.MaxStallGenerations = 30;
options.FunctionTolerance = 1e-6;
options.ConstraintTolerance = 1e-7;
options.MaxStallTime = 30;
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
ToF_min = T_min(2)-T_min(1); %days
date_start = mjd20002date(T_min(1));
date_arrival= mjd20002date(T_min(2));

fprintf ('\n start date from mercury is [%g %g %g %g %g %g] .\n',...
    date_start)
fprintf ('\n arrival date on mars is [%g %g %g %g %g %g] .\n',...
    date_arrival)

fprintf ('\n the minimum dv is [%g] km/s.\n',...
    min_dv)

t1_sec=T_min(1)*24*3600;
t2_sec=T_min(2)*24*3600;

%% plot arcs

[kep_ME,ksun] = uplanet(T_min(1), 1);
[rm_ME, v1] = sv_from_coe(kep_ME,ksun);


[kep_MA,ksun] = uplanet(T_min(2), 4);
[rm_MA, v2] = sv_from_coe(kep_MA,ksun);

ToF12_lamb_min=t2_sec-t1_sec;
[a12_min,P,E,ERROR,VI12_min,VF12_min,TPAR,THETA] = lambertMR(rm_ME,rm_MA,ToF12_lamb_min,ksun,0,0,0,0);

% FULL ORBIT IN THE YEAR OF THE ACTUAL ENCOUNTER
mjd2000_start_ME = date2mjd2000(date_start);
% one year of Mercury
oneyearME=88;
%PLOT MARS FULL ORBIT 
mjd2000_on_MA = date2mjd2000(date_arrival);
% one year of Mars
oneyearMA=687;

% for full orbit propagation
n=100;

T_MARS=linspace(mjd2000_on_MA,mjd2000_on_MA+oneyearMA,n);

T_MERC=linspace(mjd2000_start_ME,mjd2000_start_ME+oneyearME,n);

% INITIAL AND FINAL ORBIT DEFINITION
% ID NEPTUNE 8, ID MARS 4, ID MERCURY 1
% [kep,ksun] = uplanet(mjd2000, ibody)

% COMPLETE ORBITS
for k=1:n
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
plot3(rm_MA(1)./AU,rm_MA(2)./AU,rm_MA(3)./AU,'o','Color','b');
axis equal
grid minor
%legend ('full neptune','full mars','full mercury','trasnfer ME-MA','transfer MA-NE','sun','Location','southeast')
xlabel('AU')
ylabel('AU')
zlabel('AU')