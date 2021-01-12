function [fitresult, gof] = time_constant_fit(X, Y, showplot)
%time_constant_fit(X,Y)
%  Create a fit.
%  INPUT:
%      X Input : X
%      Y Output: Y
%
%  OUTPUT:
%      fitresult : a fit object representing the fit
%      gof : structure with goodness-of fit info
%

if nargin < 3
    showplot = false;
end

%% Fit:
[xData, yData] = prepareCurveData( X, Y );

% Set up fittype and options.
ft = fittype( '360*b*x', 'independent', 'x', 'dependent', 'y');
opts = fitoptions( 'Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.StartPoint = 0;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if showplot
    cla ; hold on
    title(['R^{2} = ' num2str(gof.rsquare)])
    plot(xData, yData, 'k.', 'MarkerSize', 15)
    plot(fitresult, xData, yData)
end

end