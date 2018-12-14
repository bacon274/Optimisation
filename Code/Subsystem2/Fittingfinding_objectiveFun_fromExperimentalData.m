%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting a model to experimental measured data from PVsol and MATLAB Simulink simulations %
% fitting_model to linear polynomial curve fitting regression                              %   
% looking at the goodness of fit                                                           %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num=xlsread('position_data.xlsx'); % experiment data from simulations recorded
i_0 = num(:,8);                        % (measured test simulation data & initals)
ey_0 = num(:,22);                      % (measured simualtion data test & initals)
r_0 = num(:,17);                       % (measured simulation data test & initals)

fitting_model = createFit(r_0, i_0, ey_0);
function [fitresult, gof] = createFit(R_opt_0, Irr_0, Energy_year_0)
%CREATEFIT(R_OPT_0,IRR_0,ENERGY_YEAR_0)
%  Create a fit.
%
%  Data for 'Module_Energy_Generation_fit' fit:
%      X Input : R_opt_0
%      Y Input : Irr_0
%      Z Output: Energy_year_0
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.


%% Fit: 'Module_Energy_Generation_fit'.
[xData, yData, zData] = prepareSurfaceData( R_opt_0, Irr_0, Energy_year_0 );

% Set up fittype and options.
ft = fittype( 'poly22' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Normalize = 'on';
opts.Robust = 'Bisquare';

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
figure( 'Name', 'Module_Energy_Generation_fit' );
h = plot( fitresult, [xData, yData], zData );
legend( h, 'Module_Energy_Generation_fit', 'Energy_year_0 vs. R_opt_0, Irr_0', 'Location', 'NorthEast' );
% Label axes
xlabel R_opt_0
ylabel Irr_0
zlabel Energy_year_0
grid on
view( -79.1, 12.7 );
gof
fitresult
end

