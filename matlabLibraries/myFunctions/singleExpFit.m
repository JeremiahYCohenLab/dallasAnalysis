function [fitresult, gof] = singleExpFit(coefVals, xVals, plotFlag)
%CREATEFIT(COEFVALS)
%  Create a fit.
%
%  Data for 'singleExp' fit:
%      Y Output: coefVals
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 05-Oct-2016 10:03:20

if nargin < 3
    plotFlag = 0;
end

if nargin < 2
    xVals = [];
end

%% Fit: 'singleExp'.
[xData, yData] = prepareCurveData( xVals, coefVals );

% Set up fittype and options.
ft = fittype( 'a*exp(-(1/b)*(x))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.421761282626275 0.915735525189067];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
if plotFlag == 1
    figure( 'Name', 'singleExp' );
    h = plot( fitresult, xData, yData );
    legend( h, 'coefVals', 'singleExp', 'Location', 'NorthEast' );
    % Label axes
    ylabel coefVals
    grid on
end

