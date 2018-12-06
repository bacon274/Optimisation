%Genetic Algorithm Implementation for Subsystem 1

%Genetic Algorithm Implementation for Subsystem 1

%% Global Fixed Parameters 
global a1 a2 a3 a4 a5 i1 i2 i3 i4 i5 
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

% Panel Information PS-M-NX
P1_area = 1.4;
P1_power_ub = 0.104; % note this is power /m^2 under test conditions
P1_power_lb = 0;
P1_to_g_poly_coef = [0, 3.8395e-05,   -0.0136,    1.0007] % this is the list of polynomial coefficients that map Power onto transparency

P1_cost = 600; % cost per panel

P1 = [P1_area,P1_power_ub,P1_power_lb,P1_to_g_poly_coef,P1_cost];

% Panel Information PS-CT
P2_area = 0.72;
P2_power_ub = 0.72; % note this is power /m^2 under test conditions
P2_power_lb = 0;
P2_to_g_poly_coef = [0,  0,  -0.0125, 1.0000]; % this is the list of polynomial coefficients that map Power onto CHANGE THIS FOR % TRANSPARENCY
P2_cost = 400; % cost per panel

P2 = [P2_area,P2_power_ub,P2_power_lb,P2_to_g_poly_coef,P2_cost];

%% Functions

% add linear inequality constraints

% we know from polyfit g =  -4.4103e-07*p^3 + 9.4322e-05*p^2 -0.0072*p^3 + 0.5518
% however this might need to be linear in which case g =  -0.0031*p +  0.5755

run_alg(P2)
function run_alg(P)
    global a1 a2 a3 a4 a5 i1 i2 i3 i4 i5 
    % Load Panel Properties
    P_area = P(1);
    P_power_ub = [P(2);P(2);P(2);P(2);P(2);];
    P_power_lb = [P(3);P(3);P(3);P(3);P(3);];
    P_to_g_coef = P(4:7);
    P_cost = P(8);
    
    %Objective Function
    f = @(x) objective(x,a1,a2,a3,a4,a5,i1,i2,i3,i4,i5,P_area,P_cost);
    %Constraint Fucntion
    cf = @(x) confuneq(x,a1,a2,a3,a4,a5,i1,i2,i3,i4,i5,P_to_g_coef);
    
    Number_variables = 5;
    
    [x,fval,exitflag,output] = ga(f, Number_variables,[],[],[],[],P_power_lb,P_power_ub,cf)
    
end

function z = objective(x,a1,a2,a3,a4,a5,i1,i2,i3,i4,i5,P_area,P_cost) % x1-5 = power rating  % need to add in panel information and hence work out how many panels per window
    X = [x(1) , x(2), x(3), x(4), x(5)];
    % assign costs for windows
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
    FIT = 0.386;
    %z = ((c1/(x(1)*i1*a1*FIT)) + (c2/(x(2)*i2*a2*FIT)) +(c3/(x(3)*i3*a3*FIT)) + (c4/(x(4)*i4*a4*FIT)) + (c5/(x(5)*i5*a5*FIT)))/5;
    z = (c1+c2+c3+c4+c5)/((x(1)*i1*a1*FIT)+(x(2)*i2*a2*FIT)+(x(3)*i3*a3*FIT)+(x(4)*i4*a4*FIT)+(x(5)*i5*a5*FIT));
end

function [c,ceq] = confuneq(x,a1,a2,a3,a4,a5,i1,i2,i3,i4,i5,P_to_g_coef)
    % non-linear inequality constraint 
    l = 0;
    e = 0;
    a = [a1,a2,a3,a4,a5];
    i = [i1,i2,i3,i4,i5];
    
    for j = 1:5
        g = P_to_g_coef(1)*x(j)^3 + P_to_g_coef(2)*x(j)^2 + P_to_g_coef(3)*x(j)^3 +  P_to_g_coef(4); % working out the g value from linear regression value CHANGE
        l = l + g *i(j)* a(j);% g*area*irradiance = kWh of light in  
        e = e + x(j)*i(j)*a(j); % power generated =  P * I * A
    end
    
    c1 = l - 40000000; % Light level required THIS NEEDS MORE WORK ON CONSTRAINT
    c2 = e - 500;    % Energy generated minimum 
    c = [c1;
        c2];
    ceq = [];
end