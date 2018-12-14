# Optimisation of Solar Technologies on a Domestic Property
DE4-OPTI 

This page is the github repository for the Design Engineering fourth year Optimisation Coursework. 
In this repository you will find the final report and the the code for both Subsystem 1, 2 and the whole system. 


## Prerequisites 
To run this code, you will need:

* [Matlab R2016b or above](https://uk.mathworks.com/products/matlab.html?requestedDomain=)
* [Global Optimisation Toolbox](https://uk.mathworks.com/products/global-optimization.html) 

# Subsystem 1 - Optimised configuration of Window Integrated Photovoltaics
This branch contains:
* Main optimiser code - **Subsystem_1.m**
* Exploration of the problem space - **Problem_Space_exploration.m**
* Linear regression model to determine coefficients - **Linear_regression_transparency_power.m**
* Background classes - **room_config.m , window.m**

# Subsystem 2 - Optimised Photovoltaics on roof and wall facade
This branch contains:
* Main optimiser code - **Subsystemm2.m**
* Experiment simulation models - **pv_cellModule_changing_Irr.mdl , pv_cellModule_changing_R.mdl , findMaxPower_shading.m**
* Experiment simulation test variables, result data - **position_data.xlsx**
* Variable exploration of problem space - **Exploration_of_problem_space.m**
* Curve fitting linear polynomial regression model to find objective coefficients - **Fittingfinding_objectiveFun_fromExperimentalData.m**
* Testing MultiObjective after pso - **MultiObjective_try_Subsystemm2.m**
* GA optimisation try - **run_ga_try.m**
## Authors
* **Jacob Mitchell** - *Subsystem 1* 
* **Connie Anne Dodgshon** - *Subsystem 2* 

