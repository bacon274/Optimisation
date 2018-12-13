%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Problem Space Exploration                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file creates every possible combination of window values for the
% room. 3125 combinations in total. For each combination it then calculates
% the yearly energy generated, cost, number of hours of good light let into
% the room. 
close all 
%% Compute every combination of 5 values, 0-5
N = 5;
m = dec2base(0:N^N-1, N)-'0';

% these are example combinations for testing
a = m(1,:);
b = m(30,:);

%% Create room objects with window configurations corresponding to input combination

Room_2 = room_config(b); %test examples
gval = Room_2.window_1.g;
a = m(1,:);
Room_list = [room_config(a)];
for i = 2:length(m)
    a = m(i,:);
    Room_list(i) = room_config(a);
end
%% Now need to take object and run simulations
% Simulation 1: Peak Electricity Generation

C = length(m) % number of combinations to test

%Import outdoor temperature and irridiance values (for all seasons)
Roof_I_year = 1061.6; %kWh/m^2
Wall_I_year = 744.54; % kWh/m^2

%Environment = readtable('Environmental conditions.csv'); % this is for the roof not facade
%P = zeros(N,C); % energy for each panel for each combination
%E = zeros(C,1); % energy in kWh

Cost = zeros(C,1); % cost of set up in £
Light = zeros(C,1); % total lumens into the room over the year| constraint L > 13500 Lm
Energy = zeros(C,1); %total amount of energy over the year kWh
Light_2 = zeros(C,1);
Cost_2 = zeros(C,1);
Window_list_power = zeros(C,5);
% for each combination working out the light into the room, energy
% generated, and the cost 

for i = 1:C % :length(m)
    Light_room = 0;
    Energy_room = 0;
    Room_cost = 0;
    Room = Room_list(i);
    Window_list =  [Room.window_1, Room.window_2, Room.window_3, Room.window_4, Room.window_5];
    
    %I = Environment{:,5}; % solar irridance in summer
    
    for k = 1:N
        if Window_list(k).cost ~= 0
            Room_cost = Room_cost+ Window_list(k).cost+1000;
        end

        if k > 2
            I = Wall_I_year;
        else
            I = Roof_I_year;
        end
        
        Light_panel = Window_list(k).g *I*Window_list(k).area * 105 * 683; % I THINK THIS IS INCORRECT this is light value in Lumens making it into the room over a year
        Energy_panel = (Window_list(k).power/1000)*I*Window_list(k).area;
        
        Energy_room = Energy_room + Energy_panel;
        Light_room = Light_room + Light_panel;
    end
    Window_list_power(i,:) = [Window_list.power]/1000;
    Light_2(i) = Array_Light([Window_list.power]/1000);
    Energy(i) = Energy_room;
    Light(i) = Light_room;
    Cost(i) = Room_cost;
    Cost_2(i) = Array_Cost([Window_list.power]/1000);
end

% calculate energy yeilds in kWh
%{
for j = 1:C
    for i=1:N
        hourly(i) = sum(P(i,:,j));
    end
    total = sum(hourly)/1000;
    E(j) = sum(total); % this is the energy generated in kWh that day with that combination of windows
end
%}
FIT = Energy*0.386; %this is £ generated 
Years = Cost./(FIT); 
%E_yearly = E*365;

numbers = 1:C;

T = table(numbers',Energy,FIT,Light,Years,Cost);

% inputs for linear regression
w_1 = zeros(C,1);
w_2 = zeros(C,1);
w_3 = zeros(C,1);
w_4 = zeros(C,1);
w_5 = zeros(C,1);

Array_Light([Window_list_power(1,:)]);

Array_Light([0.104,0.104,0.104,0.104,0.104]);


%% adding constraints 

g1 = 10000*10*365; %light 
g1_max = max(Light);
g2 = 500; %energy 
g2_max = max(Energy);

Gx = [g1_max g1 g1 g1_max];
Gy = [g2_max g2_max g2 g2];
Gz = [0 0 0 0];

%% Plotting

figure()

scatter3(Light_2,Energy,Years)
xlabel('Quality_light (Hours)');
ylabel('Energy Generated (kWh/year)');
zlabel('Time to pay back (years)');
set(gcf,'color','w');
hold on
title('Range space - Light vs Energy vs Time until ROI')
%patch(Gx,Gy,Gz,'green')   

rows_1 = T.Energy>g2;

vars = {'Var1','Energy','Light','Years','Cost'};
Solution_table = T(rows_1,vars);
rows_2 = Solution_table.Light>g1;
Solution_table = Solution_table(rows_2,vars)

figure()
scatter(Cost_2 , Years)
xlabel('Cost (£)');
ylabel('Time to pay back (years)');
set(gcf,'color','w');
title('Range space - Cost vs Time until ROI')

figure()
scatter(Energy , Years)
xlabel('Energy Generated (kWh/year)');
ylabel('Time to pay back (years)');
set(gcf,'color','w');
title('Range space - Energy vs Time until ROI')

figure()
scatter(Light_2 , Years)
xlabel('Quality Light (Hours)');
ylabel('Time to pay back (years)');
set(gcf,'color','w');
title('Range space - Energy vs Time until ROI')

%% Functions
function hrs = Array_Light(x)
    x;
    P_to_t_coef = [0, 38.395,   -13.6,    1.0007];
    a1= 1.44; 
    a2= 1.44; 
    a3= 0.5; 
    a4= 0.5; 
    a5= 1.8;
    % Initialise values
    Hrs_qualify = 0; % the number of hours in the day that are above the threshold amount
    
    % Parameters
    a = [a1,a2,a3,a4,a5]; % Area values for each window 
    % Yearly Irridiance values per window should these go at beginning and
    % made global?
    d = [0.4,1.5,1.5,2,2]; % distances from each window to workspace 
    I_roof = [0,0,0,0,0,0,0,0,4,107,98,84,78,67,51,24,1,0,0,0,0,0,0,0]; % irridiance values for day in january 
    I_wall = [0,0,0,0,5,80,80,100,100,80,40,30,10,7,4,2,0,0,0,0,0,0,0,0]; % irridiance values for day in january 
    
    
    % Calculate Light Levels for each hour of the day, for each window.
    % Tolling up the number of hours that qualify
    for k = 1:24
        Lux_total = 0;
        for j = 1:5
            transparency = round(polyval(P_to_t_coef,x(j)),2); % working out the g value from linear regression value CHANGE
            if j == 1 || j ==2
                I = I_roof;
            else
                I = I_wall;
            end
            
            lm = I(k)*transparency*a(j)*683; % light in lumens for window 
            Lux = lm/(9*d(j)^2); % the lux values on the table should this take area into account?
            Lux_total = Lux_total + Lux;
            
        end
        if Lux_total > 5000
            Hrs_qualify = Hrs_qualify + 1;
        end 
    end
    hrs = Hrs_qualify;
end
function c = Array_Cost(x)
    a1= 1.44; 
    a2= 1.44; 
    a3= 0.5; 
    a4= 0.5; 
    a5= 1.8;
    P_area = 1.4;
    P_cost = 400;
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
    c = c1+c2+c3+c4+c5;
end 
