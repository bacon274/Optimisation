%% Combining

% Jacob Solar window
% 1) applies initial constraints to determine all combinations of windows
% and their associated optimal energy output
% 2) All of the energy values for the different combinations are fed into
% the inner loop

%% inner loop
% 1) Energy for subsystem 2 (solar modules) to reach is = total energy
% demmand - energy output of window(for i possible combinations)
% 2) for every combination of windows connie will provide an energy output
% and total surface area used
global I
step1 = @systemPartA;
demmand = 4800;
energies = [500 400 200 300 100 600 800 700 1000];
target = step1(energies,demmand);
count0 = 0;
counts = [];
reach_energy = [15 50 100 150 180];

% finding the number of modules you need to add to the final collection of
% modules for each target energy to reach (ie for each window combination)
for i=1:length(target)
    target_energy = target(i);
end

for r=reach_energy
    if r<=target_energy
        count=count0 + 1;
        count0=count;
    end
end

for t=1:length(target)
    counts(t) = count0;
end

for c=1:(counts)
    counter=counts(c);
    final_mods_for_that_count = I(1:counter);
end

%final_modules_for_each_window_combination = ;
%final_modules_chosen = I(1:count);

%% outer loop

% target payback time defined
% difference between target and output defined (target - output)

% 1) for every output from Jacob, the inner loop with provide an output
% from Connie
% 2) Because every output from Jacob gives a daylight value and...
% 3) Because every ouput from Connie is paired to each of Jacob's combinations
% and there is a surface area and total energy  used output...
% 4) We can get energy vs daylight vs surface area
% 5) Then we can conduct multiobjective ga to trade-off daylight and
% surface area
% 6) the multiobjective is constrained to give a priority to daylight via a
% weight value the overall energy of the system will also act as a
% constraint as demmand must be met

%Pareto Front values
% 7) The output of the multiobjetcive ga will mean we can get the pareto
% front values for the correct balance of daylight to surface area
% 8) The output of the objective values correspond to the optimal
% combination of solar windows and modules
% 9) This output will have total generated energy, initial costs, payback
% time which can be compared to a target for an Analytical Target Cascading
% overall system structure approach 
% output payback time --> back to the top of the loop


% ie try a new pareto front value 