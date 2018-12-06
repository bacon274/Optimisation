%% genetic algorithm and constraints optimisation method --> running
tic 
%rng default % so that you can have reproducability
num=xlsread('position_data.xlsx');
%% step 1
ObjectiveFunction = @simple_fitness;
nvars = 2;    % Number of variables
LB=[0.3 800];         %lower bounds of variables 
UB=[0.8 2800];      %upper bounds of variables
ConstraintFunction = @simple_constraint;
[x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB, ...
    ConstraintFunction)
%% step 2
options = optimoptions(@ga,'MutationFcn',@mutationadaptfeasible);
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
%X0 = [0.5 0.5]; % Start point (row vector)
options.InitialPopulationMatrix = X0;
% Next we run the GA solver.
[x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB, ...
    ConstraintFunction,options)
