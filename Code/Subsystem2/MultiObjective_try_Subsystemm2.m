%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Multi Objective Test of Subsystem 2 - Solar Modules             %             
%    This file is a duplication of the Subsystemm2.m file, but there is the   %
%    simple testing of multiobjective ga methods to look at how the remaing   %     
%    best chosen modules could be traded off by the energy generation vs time %
%    to payback objectives; these objectives are in terms of minimizing       %    
%    the total surface area covered for a particular module layout            %    
%                            Connie Anne Dodgshon                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear variables 
clf 
close all

%% Subsystemm 2
tic 
clc 
clear all 
close all 
rng default % so that you can have reproducability
num=xlsread('position_data.xlsx');

LB=[0.055 800];         % lower bounds of variables [resistance irradiance] 
UB=[0.650 2800];        % upper bounds of variables [resistance irradiance] 

% pso parameters values 
m=2;            % number of variables 
n=37;           % population size 
wmax=0.9;       % inertia weight 
wmin=0.4;       % inertia weight 
c1=2;           % acceleration factor 

% pso main program and intial values cleaned up
maxite=5;       % set maximum number of iteration
maxrun=10;      % set maximum number of runs need to be
module_number = num(:,1);
o_fact = num(:,5);
Irr = num(:,8);     % Irr from PV sol test data set for intial elements
R_opt = num(:,17);  % R_opt from matlab test data set for intial elements

% constraint applied within the objective function for efficiency
% constraint applied below; before PSO to eliminate modules with positions within an obstruction zone on the roof or wall 

% finding the list of modules that are in obstruction zones and creating a
% list of the modules not within these areas. These remailing modules have
% module numbers and corresponding intial irradiance and initial resistance
list = [];
for i=1:length(module_number)
    if o_fact(i)<10
        list = [list; module_number(i)];
    end
end
list;

initial_irr=[];
for l=1:length(list)
    choosy=list(l);
for i=1:length(module_number)
    if module_number(i)==list(l)
        initial_irr=[initial_irr; Irr(i)];
    end
end
end
initial_irr;

initial_R_opt=[];
for l=1:length(list)
    choosy=list(l);
for i=1:length(module_number)
    if module_number(i)==list(l)
        initial_R_opt=[initial_R_opt; R_opt(i)];
    end
end
end
initial_R_opt;

%% PSO running
for run=1:maxrun     
    run

    % Energy MPP comes from MY MODEL OBJECTIVE FUNCTION so we can use the
    % swarm elements to represent a module with 52 possible places to put a
    % module - 2704 possible combinations of modules PSO needs to find the
    % best ones.
    x0 = [initial_R_opt initial_irr]; % initial values for modules from PV sol and MATLAB black box simulink test data set
    x=x0;       % initial population     
    v=0.1*x0;   % initial velocity 

    for i=1:n         
        f0(i,1)=ofun(x0(i,:));     
    end

    [fmin0,index0]=min(f0);         
    pbest=x0;               % initial pbest     
    gbest=x0(index0,:);     % initial gbest 

    %pso method
    ite=1;         
    tolerance=1; 


    while ite<=maxite && tolerance>10^-12    

        w=wmax-(wmax-wmin)*ite/maxite; % updating inertial weight

        %pso velcity updating
        for i=1:n             
            for j=1:m
                   v(i,j)=w*v(i,j)+c1*rand()*(pbest(i,j)-x(i,j)); 
            end
        end

        % pso position updating        
        for i=1:n             
            for j=1:m                 
                x(i,j)=x(i,j)+v(i,j);             
            end
        end

         % boundary violations         
         for i=1:n             
             for j=1:m                 
                 if x(i,j)<LB(j)                     
                     x(i,j)=LB(j);                 
                 elseif x(i,j)>UB(j)                     
                     x(i,j)=UB(j);                 
                 end
             end
         end

         % fitness evaluation        
         for i=1:n             
             f(i,1)=ofun(x(i,:));         
         end
         % updating pbest and fitness         
         for i=1:n             
             if f(i,1)<f0(i,1)
                 pbest(i,:)=x(i,:);                 
                 f0(i,1)=f(i,1);             
             end
         end

         [fmin,index]=min(f0);   % finding out what is the best particle         
         ffmin(ite,run)=fmin;    % remembering the best fitness         
         ffite(run)=ite;         % remembering the iteration count 

         % updating gbest and best for the fitness exploration         
         if fmin<fmin0             
             gbest=pbest(index,:);             
             fmin0=fmin;         
         end
         % looking at the tolerance         
         if ite>100;             
             tolerance=abs(ffmin(ite-100,run)-fmin0);         
         end


         % displaying iterative results         
         if ite==1             
             disp(sprintf('Iteration    Best particle    Objective fun'));         
         end
         disp(sprintf('%8g  %8g          %8.4f',ite,index,fmin0));             
         ite=ite+1; 
    end
    % new partle representing module     
    gbest;     
    fvalue=((-124.7) + (22.91*gbest(1)) + (0.03685*gbest(2))+ (1.043*gbest(1)^2) + (4.239*gbest(1)*gbest(2)) + (0.4944*gbest(2)^2))/20689;     
    fff(run)=fvalue;     
    rgbest(run,:)=gbest;     
    disp(sprintf('--------------------------------------'));
end

% pso main program 
disp(sprintf('\n')); 
disp(sprintf('*********************************************************')); 
disp(sprintf('Final Results-----------------------------')); 
[bestfun,bestrun]=min(fff) 
best_variables=rgbest(bestrun,:) 
disp(sprintf('*********************************************************')); 
time = toc ;  
%% outputs from PSO
% PSO convergence behaviour 
figure(1);
plot(ffmin(1:ffite(bestrun),bestrun),'-k'); 
xlabel('Iteration'); 
ylabel('Fitness function value'); 
title('PSO convergence behaviour');


%% MAIN SUBSYSTEM calling
Consumer_demmand = 3000; % in kWh
window_energy=158;       % in kWh, energy generated by subsystem 1 to be inputted
[Best_Modules,con] = final_mods(Consumer_demmand, window_energy, f0);
number_of_final_modules = length(Best_Modules);
     
% digitalising to make a grid and apply final inequality obstruction
% Lia = ismember(list,Best_Modules);

%% costs
% 38.3 pence per kWh is Average uk cost saving (https://www.ukpower.co.uk/home_energy/tariffs-per-unit-kwh)
per_kWh_cs =38.3;     % in pence
kWh_total = con;
c_mod_up = 900;       % �9 in pence - per module upfront cost
builder_cost = 31250; % in pence - �312.5 per 1m^2 to install (https://www.greenmatch.co.uk/blog/2014/08/what-is-the-installation-cost-for-solar-panels)

Area_mod = 0.921;     % module surface area = 0.92112m^2
Area_total = number_of_final_modules*Area_mod;

ci = (number_of_final_modules*c_mod_up) + (Area_total*builder_cost); % total upfront cost of chosen modules
yearly_saving =kWh_total*per_kWh_cs;
cmain= ci/100; % ci in �

% payback_time = installation/yearly_saving;
payback_time = ci/yearly_saving;

%% Solutions
global solution_modules
solution_modules = table(string('Particle_Swarm'),time,number_of_final_modules,payback_time,cmain,con,yearly_saving,Area_total,'VariableNames',{'Solver','Time_to_solve','Number_of_modules','Years_until_ROI','Upfront_cost','Energy_Generated','Annual_Payback','Surface_area'})


%% Multi objective on remaining modules to see impact of surface area on the model
% Aim: to see if we can use optimisation to minimise the total surface area
% of the best remaining modules that come out of the PSO algorithm and
% maximise the total energy output of the modules whilst trading off the
% payback time

% simple multi test no main constraint yet
FitnessFunction = @simple_multiobjective;
numberOfVariables = 1;
[b,mval] = gamultiobj(FitnessFunction,numberOfVariables);

size(b)
size(mval)

L = []; q = [];
Leq = []; qeq = [];
lb = 0;
ub = 500;
mul = gamultiobj(FitnessFunction,numberOfVariables,L,q,Leq,qeq,lb,ub);

options = optimoptions(@gamultiobj,'PlotFcn',{@gaplotpareto,@gaplotscorediversity});
gamultiobj(FitnessFunction,numberOfVariables,[],[],[],[],lb,ub,options);

%% functions
function f=ofun(x)

of=((124.7) + (-22.91*x(1)) + (-0.03685*x(2))+ (-1.043*x(1)^2) + (-4.239*x(1)*x(2)) + (-0.4944*x(2)^2))/20689;


c0=[]; 
c0(1)= ((38.8)/(0.55146*x(1)*x(2))*100);     % <= type constraint for module efficiency percentage

for i=1:length(c0)
    if c0(i)<6 % penalty will be 1000 if the efficiency percentage of the module is less than 6%
        c(i)=1;
    else
        c(i)=0;
    end
    
end

penalty=10000;           % penalty on each constraint violation 
f=of+penalty*sum(c); 
end

function [final_modules_chosen, energy_gen] = final_mods(Consumer_demmand, window_energy, f0)

% modules to choose - sorted descending by energy output
A = (f0*-1)-28; % *-1 t get enegy_yearly in positive from the minimization problem to maxi problem and -28 for inverter battery storage backup power (modules power themselves independantely)
[B,I] = sort(A, 'descend');
reach_energy = cumsum(B);
target = Consumer_demmand - window_energy;
count0 = 0;
for i=1:length(reach_energy)
    mod=reach_energy(i);
    if mod<=target
        count=count0 + 1;
        count0=count;
    end
end

final_modules_chosen = I(1:count);
energy_gen=reach_energy(count);
end
function y = simple_multiobjective(b)
% b is total surface area
y(1) = 27.0370.*((b./0.9210)); % objective 1 = total MPP energy from the best of the modules
y(2) = (b.*312.500*0.9210) + ((9.*b)./((b.^2).*3.46803599)); % objective 2 = payback time
end