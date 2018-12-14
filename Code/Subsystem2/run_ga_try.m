%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is the trying out of the genetic algorithm method for subsystem  % 
% 2 objective function.                                                      % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% genetic algorithm and constraints optimisation method --> running
tic 
%rng default % so that you can have reproducability
num=xlsread('position_data.xlsx');
%% step 1
ObjectiveFunction = @ofun;
nvars = 2;              % Number of variables
LB=[0.050 800];         %lower bounds of variables 
UB=[0.650 2800];        %upper bounds of variables
ConstraintFunction = @simple_constraint;
[x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB, ConstraintFunction)


%% step 1.1
ObjectiveFunction = @ofun;
nvars = 2;              % Number of variables
LB=[0.050 800];         %lower bounds of variables 
UB=[0.650 2800];        %upper bounds of variables


[x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB)



%% step 2
ObjectiveFunction = @ofun;
nvars = 2;              % Number of variables
LB=[0.050 800];         %lower bounds of variables 
UB=[0.650 2800];        %upper bounds of variables
ConstraintFunction = @simple_constraint;
options = optimoptions(@ga,'MutationFcn',@mutationadaptfeasible, 'MaxGeneration', 50, 'PopulationSize', 200);
% Next we run the GA solver.
[x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB, ...
    ConstraintFunction,options)

%% step 3 visualising
options = optimoptions(options,'PlotFcn',{@gaplotbestf,@gaplotmaxconstr}, ...
    'Display','iter');
% Next we run the GA solver.
[x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB, ...
    ConstraintFunction,options)
%% step 4 visualising with a start point
Irr = num(:,8);
R_opt = num(:,17);
X0 = [R_opt Irr];       % initial population
%X0 = [0.5 0.5];        % Start point (row vector)
options.InitialPopulationMatrix = X0;
% Next we run the GA solver.
[x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB, ...
    ConstraintFunction,options)
toc
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
% genetic algorithm and constraints optimisation method --> making
% constraints
  
function [c, ceq] = simple_constraint(x)

   %eff_constraint_percentage= ((38.8)/(0.55146*x(1)*x(2))*100);     % <= type constraint for module efficiency percentage
   %eff_constraint= ((38.8)/(0.55146*x(1)*x(2)));
   %c = @(x)(simple_fitness(x)) - (eff_constraint*simple_fitness(x));
   %c = (((124.7) + (-22.91*x(1)) + (-0.03685*x(2))+ (-1.043*x(1)^2) + (-4.239*x(1)*x(2)) + (-0.4944*x(2)^2))/20689) - (eff_constraint*(((124.7) + (-22.91*x(1)) + (-0.03685*x(2))+ (-1.043*x(1)^2) + (-4.239*x(1)*x(2)) + (-0.4944*x(2)^2))/20689));

   c= ((38.8)/(0.55146*x(1)*x(2)));

   ceq = [];
end