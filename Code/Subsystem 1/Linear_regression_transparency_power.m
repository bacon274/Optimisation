%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Linear Regression                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is run to find the equation that maps the transparency onto the
% power variables. This is so it can be modelled as a continuous variable
% in the solvers in Subsystem_1.m

% Power and transparency values
p = [0.104,0.088,0.062,0.042,0]
t = [0,0.1,0.3,0.5,1];

p2 = [0.072,0.064,0.056,0.048,0.040,0]; 
t2 = [0.1,0.2,0.3,0.4,0.50,1];  

% function fitting and plot results to check fit
z = fit_function(p,t)

transparency = polyval(z,0.042)

function z = fit_function(p,g)
    for i = 1:2

        z = polyfit(p,g,i)
        p1 = linspace(0,p(1));
        g1 = polyval(z,p1);
        
        figure()
        scatter(p,g)
        hold on
        plot(p1,g1)
        xlabel('Power (kWh/m^2)')
        ylabel('Transparency (%)')
        title('Linear Regression of Transparency and Power')
    end
end 



 