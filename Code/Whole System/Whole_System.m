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

