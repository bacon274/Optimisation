%Test of variable selection

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
P1_power_ub = 104; % note this is power /m^2 under test conditions
P1_power_lb = 0;
P1_to_g_poly_coef = [38.3949  -13.6191    1.0007]; % this is the list of polynomial coefficients that map Power onto transparency

P1_cost = 600; % cost per panel

P1 = [P1_area,P1_power_ub,P1_power_lb,P1_to_g_poly_coef,P1_cost];

% Panel Information PS-CT
P2_area = 0.72;
P2_power_ub = 720; % note this is power /m^2 under test conditions
P2_power_lb = 0;
P2_to_g_poly_coef = [0,  0,  -12.5, 1.0000]; % this is the list of polynomial coefficients that map Power onto CHANGE THIS FOR % TRANSPARENCY
P2_cost = 400; % cost per panel

P2 = [P2_area,P2_power_ub,P2_power_lb,P2_to_g_poly_coef,P2_cost];

%% Functions
x = [ 0    0.1040         0.1040         0         0];
P_area = P1_area; 
P_cost = P1_cost; 
P_to_t_coef = P1_to_g_poly_coef;

f = objective(x,a1,a2,a3,a4,a5,i1,i2,i3,i4,i5,P_area,P_cost);
    %Constraint Fucntion
cf = confuneq(x,a1,a2,a3,a4,a5,i1,i2,i3,i4,i5,P_to_t_coef);

function Years = objective(x,a1,a2,a3,a4,a5,i1,i2,i3,i4,i5,P_area,P_cost) % x1-5 = power rating  % need to add in panel information and hence work out how many panels per window
    X = [x(1) , x(2), x(3), x(4), x(5)]
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
    Years = (c1+c2+c3+c4+c5)/((x(1)*i1*a1*FIT)+(x(2)*i2*a2*FIT)+(x(3)*i3*a3*FIT)+(x(4)*i4*a4*FIT)+(x(5)*i5*a5*FIT))
end

function [c,ceq] = confuneq(x,a1,a2,a3,a4,a5,i1,i2,i3,i4,i5,P_to_t_coef)
    % non-linear inequality constraint 
    
    % Initialise values
    Hrs_qualify = 0; % the number of hours in the day that are above the threshold amount
    e = 0;
    
    % Parameters
    a = [a1,a2,a3,a4,a5]; % Area values for each window 
    i = [i1,i2,i3,i4,i5]; % Yearly Irridiance values per window
    d = [0.4,1.5,1.5,2,2]; % distances from each window to workspace 
    I_roof = [0,0,0,0,0,0,0,0,4,107,98,84,78,67,51,24,1,0,0,0,0,0,0,0]; % irridiance values for day in january 
    I_wall = [0,0,0,0,5,80,80,100,100,80,40,30,10,7,4,2,0,0,0,0,0,0,0,0]; % irridiance values for day in january 
    
    
    % Calculate Light Levels for each hour of the day, for each window.
    % Tolling up the number of hours that qualify
    for k = 1:24
        Lux_total = 0;
        for j = 1:5
            transparency = round(polyval(P_to_t_coef,x(j)),2)
           % transparency = P_to_t_coef(1)*x(j)^3 + P_to_t_coef(2)*x(j)^2 + P_to_t_coef(3)*x(j)^3 +  P_to_t_coef(4) % working out the g value from linear regression value CHANGE
            if j == 1 || j ==2
                I = I_roof;
            else
                I = I_wall;
            end
            
            lm = I(k)*transparency*a(j)*683; % light in lumens for window (
            Lux = lm/(9*d(j)^2); % the lux values on the table should this take area into account?
            Lux_total = Lux_total + Lux;
        end
        
        if Lux_total > 10000
            Hrs_qualify = Hrs_qualify + 1;
        end 
    end
    
    % Calculate the Yearly Energy Production Values
    for j = 1:5
        e = e + x(j)*(i(j))*a(j); % power generated =  P (kW) * I(kwh/m^2) * A (m^2)?
    end
    
    Hrs_qualify
    e
    
    c1 = -Hrs_qualify + 5;  % 4 hours of good natural light per day 
    c2 = 1000 -e ;          % Energy generated minimum 
    c = [c1;
        c2];
    ceq = [];
end