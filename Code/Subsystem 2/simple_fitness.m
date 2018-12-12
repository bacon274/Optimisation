% genetic algorithm and constraints optimisation method --> making function
function y = simple_fitness(x)
    y =((124.7) + (-22.91*x(1)) + (-0.03685*x(2))+ (-1.043*x(1)^2) + (-4.239*x(1)*x(2)) + (-0.4944*x(2)^2))/20689;
    

