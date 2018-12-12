tic 
clc 
clear all 
close all 
rng default % so that you can have reproducability
num=xlsread('position_data.xlsx');

LB=[0.055 800];         %lower bounds of variables 
UB=[0.650 2800];      %upper bounds of variables

% pso parameters values 
m=2;            % number of variables 
n=52;          % population size 
wmax=0.9;       % inertia weight 
wmin=0.4;       % inertia weight 
c1=2;           % acceleration factor 
%c2=2;

% pso main program
maxite=5;    % set maximum number of iteration
maxrun=10;      % set maximum number of runs need to be

for run=1:maxrun     
    run

    % pso initialisation
    %for i=1:n         
     %   for j=1:m
      %      x0(i,j)=round(LB(j)+rand()*(UB(j)-LB(j)));         
       % end 
    %end
    Irr = num(:,8);
    R_opt = num(:,17);
    x0 = [R_opt Irr];
    x=x0;       % initial population     
    v=0.1*x0;   % initial velocity 

    for i=1:n         
        f0(i,1)=ofun(x0(i,:));     
    end

    [fmin0,index0]=min(f0);         
    pbest=x0;               % initial pbest     
    gbest=x0(index0,:);     % initial gbest 

    %pso algorithim
    ite=1;         
    tolerance=1; 


    while ite<=maxite && tolerance>10^-12    

        w=wmax-(wmax-wmin)*ite/maxite; % update inertial weight

        %pso vvelcity updates
        for i=1:n             
            for j=1:m
                   v(i,j)=w*v(i,j)+c1*rand()*(pbest(i,j)-x(i,j)); 
            end
        end

        % pso position update         
        for i=1:n             
            for j=1:m                 
                x(i,j)=x(i,j)+v(i,j);             
            end
        end

         % handling boundary violations         
         for i=1:n             
             for j=1:m                 
                 if x(i,j)<LB(j)                     
                     x(i,j)=LB(j);                 
                 elseif x(i,j)>UB(j)                     
                     x(i,j)=UB(j);                 
                 end
             end
         end

         % evaluating fitness         
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

         [fmin,index]=min(f0);   % finding out the best particle         
         ffmin(ite,run)=fmin;    % storing best fitness         
         ffite(run)=ite;         % storing iteration count 

         % updating gbest and best fitness         
         if fmin<fmin0             
             gbest=pbest(index,:);             
             fmin0=fmin;         
         end
         % calculating tolerance         
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
    % pso algorithm     
    gbest;     
    fvalue=((-124.7) + (22.91*gbest(1)) + (0.03685*gbest(2))+ (1.043*gbest(1)^2) + (4.239*gbest(1)*gbest(2)) + (0.4944*gbest(2)^2))/20689;

    %fvalue=10*(gbest(1)-1)^2+20*(gbest(2)-2)^2+30*(gbest(3)-3)^2;     
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
toc   
%% outputs
% PSO convergence characteristic 
figure(1)
plot(ffmin(1:ffite(bestrun),bestrun),'-k'); 
xlabel('Iteration'); 
ylabel('Fitness function value'); 
title('PSO convergence characteristic')
%figure(2)
%scatter3(x(:,1),x(:,2), f0)

global I
% modules to choose - sorted descending by energy output
A = (f0*-1)-28; % *-1 t get enegy_yearly in positive from the minimization problem to maxi problem and -28 for inverter battery storage backup power (modules power themselves independantely)
[B,I] = sort(A, 'descend');
reach_energy = cumsum(A);
target = 3500;
count0 = 0;
for i=1:length(reach_energy)
    mod=reach_energy(i);
    if mod<=target
        count=count0 + 1;
        count0=count;
    end
end

final_modules_chosen = I(1:count);
%% costs

% 12.827 pence per kWh is Average uk cost saving (https://www.ukpower.co.uk/home_energy/tariffs-per-unit-kwh)
per_kWh_cs =38.6;%12.827; % in pence
kWh_total = reach_energy(count);
c_mod_up = 900; % £9 in pence - per module upfront cost
%c_build_sal = 1225; % £12.25  in pernce - per hour builder's average salary in uk (https://www.payscale.com/research/UK/Job=Builder/Hourly_Rate)
builder_cost = 31250; % in pence - £312.5 per 1m^2 to install (https://www.greenmatch.co.uk/blog/2014/08/what-is-the-installation-cost-for-solar-panels)

Area_mod = 0.921; % 0.92112m^2
Area_total = count*Area_mod;

ci = (count*c_mod_up) + (Area_total*builder_cost);
yearly_saving =kWh_total*per_kWh_cs;
cmain= ci/100;

%payback_time = installation/yearly_saving;
payback_time = ci/yearly_saving;



%% Multi objective on remaining modules to see impact of surface area on the model
% Aim: to see if we can use optimisation to minimise the total surface area
% of the best remaining modules that come out of the PSO algorithm and
% maximise the total energy output of the modules
% constraint is to hit the efficiency required

% Plot two objective functions on the same axis
count;
Av_mod_energy = mean(final_modules_chosen);
Av_yearly_saving = Av_mod_energy*per_kWh_cs;
module_SA = Area_mod; % 0.92112m^2;
builder_pounds = builder_cost/100; %31250 pence;
max_SA = count*module_SA;

b = 0:module_SA:max_SA;
%m1 = (Av_mod_energy.*(b./module_SA)); % total MPP energy from the best of the modules
%m2 = ((b.*builder_pounds) + 9.*(b./module_SA))./Av_yearly_saving; % Average payback time: total investment costs of modules/ Average yearly cost saving
m1 = 27.0370.*((b./0.9210)); % total MPP energy from the best of the modules
m2 = (b.*312.500*0.9210) + ((9.*b)./((b.^2).*3.46803599)); %payback time and discounting from market and builder
figure(3)
plot(b,m1);
hold on;
plot(b,m2,'r');
grid on;
title('Plot of objectives ''Energy'' and ''Cost''');
legend('Energy','Cost')
%% multi2 test
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













%% for multiobjective
%o_fact = num(:,5);
%restriction = o_fact/100;
%eff = ((38.8)/(0.55146*x(1)*x(2))*100); % eff is the general efficiency for each module based on position relative to irr and shading
%efficiency_o_fact = (eff) - (eff*restriction); % efficiency_o_fact is the overall efficiency that corrects the eff for each module due to the impact of an obstruction on the building in the proposed position of the module

% the efficiency_o_fact for each module should be maximised the o_fact
% effects the efficiency due to the obstruction having to break the circuit
% and requiring a separate battery/ circuit module configuration. The
% o_fact value comes from the energy loss due to the probability of having to use another circuit configuration instead of staying constistant with the most efficient series- parrel config with is simulatied in pv sol for all modules duplicated assuming there is no physical design obstruction on the roof/ wall position coordinate

%% display as on or off modules grid
module_number = num(:,1);
roof_wall_x = num(:,2); % centre of module x coord
roof_wall_y = num(:,3); % centre of module y coord

%[X,Y,Z] = meshgrid(roof_wall_x,roof_wall_y, module_number);
%figure(2)
%surf(roof_wall_x, roof_wall_y, module_number)