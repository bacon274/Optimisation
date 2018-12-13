%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Whole System Analysis                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file runs the whole system analysis by running subsystem 1, obtaining
% an optimised solution for the window configuration. This is shared as a 
% global variable between files. The energy required of subsystem 2 is 
% calculated and shared as a global variable.Subsystem 2 is then run and 
% it calculates its own optimal solution. 

%The outputs of both solutions are then outputted at the end of the file
%for inspection

clc 
close all 
clear all 

global solution_table target_energy solution_table_2

Subsystem_1
Total_energy = 3000; 
window_energy = solution_table{:,'Energy_Generated'}
target_energy = 3000 - window_energy;
FinalSubsystem2
solution_table_2
solution_table


Total_Energy = solution_table.Energy_Generated + solution_table_2.Energy_Generated
Total_Cost = solution_table.Upfront_cost + solution_table_2.Upfront_cost
Total_Years = Total_Cost/(Total_Energy*0.386)

