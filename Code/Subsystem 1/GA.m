%%Genetic Algorithm Implementation for Subsystem 1
%area and irridiance are fixed parameters 
a1 = 1.44; 
a2 = 1.44; 
a3= 0.15; 
a4= 0.15; 
a5 =1.05;

i1=1061.16; % this is kWh/m^2
i2=1061.16;
i3=744.54;
i4=744.54;
i5=744.54;

f = @(x) objective(x,a1,a2,a3,a4,a5,i1,i2,i3,i4,i5);

cf = @(x) confuneq(x,a1,a2,a3,a4,a5,i1,i2,i3,i4,i5);

Number_variables = 5;

ub = [0.160;0.160;0.160;0.160;0.160;]; % this is power in kilowatts
lb = [0;0;0;0;0;];

[x,fval,exitflag,output] = ga(f, Number_variables,[],[],[],[],lb,ub,cf)


% add linear inequality constraints

% we know from polyfit g =  -4.4103e-07*p^3 + 9.4322e-05*p^2 -0.0072*p^3 + 0.5518
% however this might need to be linear in which case g =  -0.0031*p +  0.5755



function z = objective(x,a1,a2,a3,a4,a5,i1,i2,i3,i4,i5) % x1-5 = power rating 
    % assign costs for windows
    if x(1) == 0 
        c1 = 0;
    else
        c1 = 600 + 500;
    end
    
    if x(2) == 0
        c2 = 0;
    else
        c2 = 600 + 500;
    end
    
    if x(3) == 0
        c3 = 0;
    else
        c3 = 600 + 500;
    end
    
    if x(4) == 0
        c4 = 0;
    else
        c4 = 600 + 500;
    end
    
    if x(5) == 0
        c5 = 0;
    else
        c5 = 600 + 500;
    end
    FIT = 0.386;
    %z = ((c1/(x(1)*i1*a1*FIT)) + (c2/(x(2)*i2*a2*FIT)) +(c3/(x(3)*i3*a3*FIT)) + (c4/(x(4)*i4*a4*FIT)) + (c5/(x(5)*i5*a5*FIT)))/5;
    z = (c1+c2+c3+c4+c5)/((x(1)*i1*a1*FIT)+(x(2)*i2*a2*FIT)+(x(3)*i3*a3*FIT)+(x(4)*i4*a4*FIT)+(x(5)*i5*a5*FIT));
end

function [c,ceq] = confuneq(x,a1,a2,a3,a4,a5,i1,i2,i3,i4,i5)
    % non-linear inequality constraint 
    l = 0;
    e = 0;
    a = [a1,a2,a3,a4,a5];
    i = [i1,i2,i3,i4,i5];
    
    for j = 1:5
        g = -4.4103e-07*x(j)^3 + 9.4322e-05*x(j)^2 -0.0072*x(j)^3 + 0.5518; % working out the g value from linear regression value
        l = l + g *i(j)* a(j);% g*area*irradiance = kWh of light in  
        e = e+ x(j)*i(j)*a(j); % power generated =  P * I * A
    end
      
    c1 = l - 40000000;
    c2 = e - 500;
    c = [c1;
        c2];
    ceq = [];
end