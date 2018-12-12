% genetic algorithm and constraints optimisation method --> making
% constraints
  
function [c, ceq] = simple_constraint(x)

   %eff_constraint_percentage= ((38.8)/(0.55146*x(1)*x(2))*100);     % <= type constraint for module efficiency percentage
   eff_constraint= ((38.8)/(0.55146*x(1)*x(2)));
   %c = @(x)(simple_fitness(x)) - (eff_constraint*simple_fitness(x));
   c = (((124.7) + (-22.91*x(1)) + (-0.03685*x(2))+ (-1.043*x(1)^2) + (-4.239*x(1)*x(2)) + (-0.4944*x(2)^2))/20689) - (eff_constraint*(((124.7) + (-22.91*x(1)) + (-0.03685*x(2))+ (-1.043*x(1)^2) + (-4.239*x(1)*x(2)) + (-0.4944*x(2)^2))/20689));

   ceq = [];