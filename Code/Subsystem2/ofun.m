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
