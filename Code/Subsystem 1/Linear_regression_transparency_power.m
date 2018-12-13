% linear regression Power and g value 
p = [0.104,0.088,0.062,0.042,0]
t = [0,0.1,0.3,0.5,1];


p2 = [0.072,0.064,0.056,0.048,0.040,0]; 
t2 = [0.1,0.2,0.3,0.4,0.50,1];
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
    end
end 



 