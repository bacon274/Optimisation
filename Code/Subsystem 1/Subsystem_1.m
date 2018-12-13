%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Subsystem 1 - BIPV Windows                      %
%   This file is the implementation of Subsystem 1 optimisation analysis  %
%      it runs independently or when called with the whole system file    %
%                              Jacob Mitchell                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear variables 
clf 
close all

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Global Fixed Parameters                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global a1 a2 a3 a4 a5 i1 i2 i3 i4 i5 Results_Array_Energy Results_Array_Years Results_Array_Hrs_light Results_Array_Cost solution solution_table

% Area of Windows (m^2)
a1= 1.44; 
a2= 1.44; 
a3= 0.5; 
a4= 0.5; 
a5= 1.8;
% Yearly Light Irridance on Windows (kWh/m^2)
i1=1061.16;  
i2=1061.16;
i3=744.54;
i4=744.54;
i5=744.54;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Results Arrays for plotting outcomes of solvers          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results_Array_Energy = [];
Results_Array_Years = [];
Results_Array_Hrs_light = [];
Results_Array_Cost = [];
solution_table = table();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Panel information                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Panel Information PS-M-NX
P1_area = 1.4;
P1_power_ub = 0.104; % note this is power /m^2 under test conditions
P1_power_lb = 0;
P1_to_g_poly_coef = [0, 38.395,   -13.6,    1.0007]; % this is the list of polynomial coefficients that map Power onto transparency
P1_cost = 400; % cost per panel
P1 = [P1_area,P1_power_ub,P1_power_lb,P1_to_g_poly_coef,P1_cost];

% Panel Information PS-CT
P2_area = 0.72;
P2_power_ub = 0.072; % note this is power /m^2 under test conditions
P2_power_lb = 0;
P2_to_g_poly_coef = [0,  0,  -12.5, 1.0000]; % this is the list of polynomial coefficients that map Power onto CHANGE THIS FOR % TRANSPARENCY
P2_cost = 400; % cost per panel
P2 = [P2_area,P2_power_ub,P2_power_lb,P2_to_g_poly_coef,P2_cost];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Solver Functions                              %
%               Uncomment the functions below to run each solver          %           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Both Panel type 1 and 2 can be tested however, initial tests found P1 to 
% be more cost effective

run_ga(P1)                         % Genetic Algorithm
%run_fmincon(P1)                   % Gradient Based Constrained Optimiser
%run_swarm(P1)                     % Particle Swarm Analysis 
%run_ga2(P1)                       % Genetic Algorithm w/constraints in obj
%run_multiobjectga(P1)             % Multi Objective Genetic Algorithm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Post Optimal Analysis                           %
%                Uncomment the functions below to run                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sensitivity([0.104;0.104;0.104	0.104;0.104] ,P1) 

%Test([0,0,0,0,0.104],P1_to_g_poly_coef,P1_area,P2_cost)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Functions                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm Functions
function run_fmincon(P)
    P_area = P(1);
    P_power_ub = [P(2);P(2);P(2);P(2);P(2);];
    P_power_lb = [0;0;0;0;0;];
    P_to_t_coef = P(4:7);
    P_cost = P(8);
    
    %Objective Function
    f = @(x) objective(x,P_area,P_cost);
    %Constraint Fucntion
    cf = @(x) confuneq(x,P_to_t_coef,P_area,P_cost);
    x0 = [0 0 0 0 0];
    tic
    [x,fval,exitflag,output] = fmincon(f,x0,[],[],[],[],P_power_lb,P_power_ub,cf)
    time = toc
    Test(x,P_to_t_coef,P_area,P_cost,time,'fmincon')
end
function run_ga(P)
    % Load Panel Properties
    P_area = P(1);
    P_power_ub = [P(2);P(2);P(2);P(2);P(2);];
    P_power_lb = [0;0;0;0;0;];
    P_to_t_coef = P(4:7);
    P_cost = P(8);
    
    %Objective Function
    f = @(x) objective(x,P_area,P_cost);
    
    %Constraint Fucntion
    cf = @(x) confuneq(x,P_to_t_coef,P_area,P_cost);
    
    Number_variables = 5;
    X0 = [0.104 0.104 0.104 0.104 0.104]; % start point
    options.InitialPopulationMatrix = X0; 
    tic 
    [x,fval,exitflag,output] = ga(f, Number_variables,[],[],[],[],P_power_lb,P_power_ub,cf,options)
    time = toc;
    if exitflag == -2
        Clean()
        run_ga(P)
    else
        Test(x,P_to_t_coef,P_area,P_cost,time,'Genetic Algorithm')
        Clean()
    end
    
    
end
function run_ga2(P)
    % Load Panel Properties
    P_area = P(1);
    P_power_ub = [P(2);P(2);P(2);P(2);P(2);];
    P_power_lb = [0;0;0;0;0;];
    P_to_t_coef = P(4:7);
    P_cost = P(8);
    
    %Objective Function
   % f = @(x) objective(x,P_area,P_cost);
    f = @(x) SO(x,P_to_t_coef,P_area,P_cost);
    
    
    %Constraint Fucntion
    %cf = @(x) confuneq(x,P_to_t_coef,P_area,P_cost);
    
    Number_variables = 5;
    X0 = [0 0 0 0 0]; % start point
    options.InitialPopulationMatrix = X0; 
    tic 
    [x,fval,exitflag,output] = ga(f, Number_variables,[],[],[],[],P_power_lb,P_power_ub,[],options)
    time = toc;
    if exitflag == -2
        Clean()
        run_ga(P)
    else
        Test(x,P_to_t_coef,P_area,P_cost,time,'Genetic Algorithm')
        Clean()
    end
    
    
end
function run_swarm(P)
    P_area = P(1);
    P_power_ub = [P(2);P(2);P(2);P(2);P(2);];
    P_power_lb = [0;0;0;0;0;];
    P_to_g_coef = P(4:7);
    P_cost = P(8);
    
    %Objective Function
    f = @(x) SO(x,P_to_g_coef,P_area,P_cost);
    
    Number_variables = 5;
    X0 = [0 0 0 0 0]; % start pointNumber_variables
    tic
    [x,fval,exitflag] = particleswarm(f,Number_variables,P_power_lb,P_power_ub)
    time = toc;
    Test(x,P_to_g_coef,P_area,P_cost,time,'Particle Swarm')
    Clean()
    
end
function run_multiobjectga(P)
    P_area = P(1);
    P_power_ub = [P(2);P(2);P(2);P(2);P(2);];
    P_power_lb = [0;0;0;0;0;];
    P_to_t_coef = P(4:7);
    P_cost = P(8);
    
    %Objective Function
    f = @(x) [Array_cost(x,P_area,P_cost); -Array_energy(x); -Light(x,P_to_t_coef)]; %objective(x,P_area,P_cost) SO(x,P_to_t_coef,P_area,P_cost)
    %Constraint Fucntion
    %cf = @(x) confuneq(x,P_to_g_coef,P_area,P_cost);
    
    Number_variables = 5;
    X0 = [0 0 0 0 0]; % start point
    %options.InitialPopulationMatrix = X0;
    options = optimoptions(@gamultiobj,'PlotFcn',{@gaplotpareto,@gaplotscorediversity});
    tic 
    [x,fval,exitflag,output] = gamultiobj(f,Number_variables,[],[],[],[],P_power_lb,P_power_ub,options)
    time = toc;
    for i = 1:length(x)
        x(i,:)
        T = Test(x(i,:),P_to_t_coef,P_area,P_cost,time,'Multi Objective Genetic Algorithm')
    end
   
end

% Objective and Constraint Functions
function [c,ceq] = confuneq(x,P_to_t_coef,P_area,P_cost)
    global Results_Array_Energy Results_Array_Hrs_light Results_Array_Cost Results_Array_Years
    % non-linear inequality constraint 
    %Light Hours of the day that qualify
    Hrs_qualify = Light(x,P_to_t_coef);
    %Energy generated by Array
    Energy_generated = Array_energy(x);
    %Array cost
    Cost = Array_cost(x,P_area,P_cost);
   
    c1 = 6 - Hrs_qualify;    % 5 hours of good natural light per day 
    c2 = 100 - Energy_generated ;  % Energy generated minimum 
    c3 = Cost-4000; % cost belw £5000
   

    c = [c1;
         c2
         c3];
    ceq = [];
    
    %Add resuls to Results Array 
    Results_Array_Energy = [Results_Array_Energy, Energy_generated];
    Results_Array_Hrs_light = [Results_Array_Hrs_light, Hrs_qualify];
    Results_Array_Cost = [Results_Array_Cost, Cost];
    Results_Array_Years = [Results_Array_Years, objective(x,P_area,P_cost)];
end
function z = objective(x,P_area,P_cost) % x1-5 = power rating  % need to add in panel information and hence work out how many panels per window
    % assign costs for windows
    Cost = Array_cost(x,P_area,P_cost);
    Energy_generated = Array_energy(x);
    Annual_payback = FIT_payback(Energy_generated);
    Years_to_payback = Years(Cost,Annual_payback);
    z = Years_to_payback;
end

%Swarm Objective function including constraints 
function z = SO(x,P_to_t_coef,P_area,P_cost)
% assign costs for windows
    global Results_Array_Energy Results_Array_Hrs_light Results_Array_Cost Results_Array_Years

    Cost = Array_cost(x,P_area,P_cost);
    Energy_generated = Array_energy(x);
    Annual_payback = FIT_payback(Energy_generated);
    Years_to_payback = Years(Cost,Annual_payback);
    Hrs_qualify = Light(x,P_to_t_coef);
    penalty = 0;
    Results_Array_Energy = [Results_Array_Energy, Energy_generated];
    Results_Array_Hrs_light = [Results_Array_Hrs_light, Hrs_qualify];
    Results_Array_Cost = [Results_Array_Cost, Cost];
    Results_Array_Years = [Results_Array_Years, Years_to_payback];
    
    
    if  Hrs_qualify < 6  ||Cost > 4000  || Energy_generated < 100
        penalty = 1000;
    end
    
    
    z = Years_to_payback + penalty;

end

% Calculation Functions:
function e = Array_energy(x)
    global a1 a2 a3 a4 a5 i1 i2 i3 i4 i5
    e = ((x(1)*i1*a1)+(x(2)*i2*a2)+(x(3)*i3*a3)+(x(4)*i4*a4)+(x(5)*i5*a5));
end
function f = FIT_payback(e)
    f = e*0.386;
end
function y = Years(c,f)
    y = c/f;
end
function c = Array_cost(x,P_area,P_cost)
    global a1 a2 a3 a4 a5
    if x(1) == 0 
        c1 = 0;
    else
        panels = ceil(a1/P_area);
        c1 = P_cost*panels + 500;
    end
    
    if x(2) == 0
        c2 = 0;
    else
        panels = ceil(a2/P_area);
        c2 = P_cost*panels + 500;
    end
  
    if x(3) == 0
        c3 = 0;
    else
        panels = ceil(a3/P_area);
        c3 = P_cost*panels + 500;
    end
    
    if x(4) == 0
        c4 = 0;
    else
        panels = ceil(a4/P_area);
        c4 = P_cost*panels + 500;
    end
    
    if x(5) == 0
        c5 = 0;
    else
        panels = ceil(a5/P_area);
        c5 = P_cost*panels + 500;
    end
    c = c1+c2+c3+c4+c5;
end 
function hrs = Light(x,P_to_t_coef)
    global a1 a2 a3 a4 a5 
    % Initialise values
    Hrs_qualify = 0; % the number of hours in the day that are above the threshold amount
    
    % Parameters
    a = [a1,a2,a3,a4,a5]; % Area values for each window 
    % Yearly Irridiance values per window should these go at beginning and
    % made global?
    d = [0.4,1.5,1.5,2,2]; % distances from each window to workspace 
    I_roof = [0,0,0,0,0,0,0,0,4,107,98,84,78,67,51,24,1,0,0,0,0,0,0,0]; % irridiance values for day in january 
    I_wall = [0,0,0,0,5,80,80,100,100,80,40,30,10,7,4,2,0,0,0,0,0,0,0,0]; % irridiance values for day in january 
    
    
    % Calculate Light Levels for each hour of the day, for each window.
    % Tolling up the number of hours that qualify
    for k = 1:24
        Lux_total = 0;
        for j = 1:5
            transparency = round(polyval(P_to_t_coef,x(j)),2); % working out the g value from linear regression value CHANGE
            if j == 1 || j ==2
                I = I_roof;
            else
                I = I_wall;
            end
            
            lm = I(k)*transparency*a(j)*683; % light in lumens for window 
            Lux = lm/(9*d(j)^2); % the lux values on the table should this take area into account?
            Lux_total = Lux_total + Lux;
            
        end
        if Lux_total > 30000
            Hrs_qualify = Hrs_qualify + 1;
        end 
    end
    hrs = Hrs_qualify;
end
function l = Light2(x,P_to_t_coef)
    global a1 a2 a3 a4 a5 
    % Initialise values
     % the number of hours in the day that are above the threshold amount
    
    % Parameters
    a = [a1,a2,a3,a4,a5]; % Area values for each window 
    % Yearly Irridiance values per window should these go at beginning and
    % made global?
    d = [0.4,1.5,1.5,2,2]; % distances from each window to workspace 
    I_roof = [0,0,0,0,0,0,0,0,4,107,98,84,78,67,51,24,1,0,0,0,0,0,0,0]; % irridiance values for day in january 
    I_wall = [0,0,0,0,5,80,80,100,100,80,40,30,10,7,4,2,0,0,0,0,0,0,0,0]; % irridiance values for day in january 
    Lux_array = [];
    
    % Calculate Light Levels for each hour of the day, for each window.
    % Tolling up the number of hours that qualify
    for k = 1:24
        Lux_total = 0;
        for j = 1:5
            transparency = round(polyval(P_to_t_coef,x(j)),2); % working out the g value from linear regression value CHANGE
            if j == 1 || j ==2
                I = I_roof;
            else
                I = I_wall;
            end
            
            lm = I(k)*transparency*a(j)*683; % light in lumens for window 
            Lux = lm/(9*d(j)^2); % the lux values on the table should this take area into account?
            Lux_total = Lux_total + Lux;
            
        end
        Lux_array = [Lux_array; Lux_total;];
    end 
    l = sum(Lux_array);

end

% Post analysis
function sensitivity(x,P)
    P_area = P(1);
    P_to_t_coef = P(4:7);
    P_cost = P(8);
    solver = 'Sensitivity Test';
    time = 1
    Test(x,P_to_t_coef,P_area,P_cost,time,solver)
    for i = 1:5
        X = x;
        X(i) = X(i)- X(i)/10;
        Test(X,P_to_t_coef,P_area,P_cost,time,solver)
    end
end

%Test Output from solver
function [T] = Test(x,P_to_t_coef,P_area,P_cost,time,solver)
    global Results_Array_Energy Results_Array_Years Results_Array_Hrs_light  solution_table
   
    Cost_pounds = Array_cost(x,P_area,P_cost);
    Energy_generated_kWh = Array_energy(x);
    Annual_payback_pounds = FIT_payback(Energy_generated_kWh);
    Years_to_payback = Years(Cost_pounds,Annual_payback_pounds);
    Hrs_light_qualify = Light(x,P_to_t_coef);
    
    T = table(string(solver),time,x(1),x(2),x(3),x(4),x(5),Years_to_payback,Cost_pounds,Energy_generated_kWh,Annual_payback_pounds,Hrs_light_qualify,'VariableNames',{'Solver','Time_to_solve', 'W1','W2','W3','W4','W5','Years_until_ROI','Upfront_cost','Energy_Generated','Annual_Payback','Hrs_of_quality_light'});
    
    solution_table = [solution_table; T;];
    GX = [[0,1,1,0];     [1.5,2.5,2.5,1.5]; [0,1,1,0];     [1.5,2.5,2.5,1.5]; [0.75,1.75,1.75,0.75];];
    GY = [[3,3,4.4,4.4]; [3,3,4.4,4.4];     [2,2,2.5,2.5]; [2,2,2.5,2.5];    [0,0,1.8,1.8];];
    figure()
    hold on
    for i = 1:5
        trans = (1-(x(i)/0.104))*100;
        patch(GX(i,:),GY(i,:),trans);
    end
    cb = colorbar;
    cb.Label.String = 'Transparency %' ;
    set(gcf,'color','w');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    
    
    figure()
    scatter3(Results_Array_Energy,Results_Array_Hrs_light,Results_Array_Years)
    xlabel('Energy (kWh/year)')
    ylabel('Work Quality Light (Hrs/day)')
    zlabel('Years to pay back')
    hold on 
    scatter3(Energy_generated_kWh, Hrs_light_qualify,Years_to_payback,[100],'red','filled')
    set(gcf,'color','w');
end

% House Keeping
function Clean()
    global Results_Array_Energy Results_Array_Years Results_Array_Hrs_light Results_Array_Cost
    Results_Array_Energy = [];
    Results_Array_Years  = [];
    Results_Array_Hrs_light  = [];
    Results_Array_Cost = [];
end
