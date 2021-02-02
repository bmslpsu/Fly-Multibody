function [fitresult, gof] = fitSS(x, y, f, n_detrend, showplot)
%fitSS(TIME,Y)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : time
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%

if nargin < 5
   showplot = true;
   if nargin < 4
       n_detrend = 0;
   end
end

%% Fit
[xData, yData] = prepareCurveData( x, y );

yData = detrend(yData, 10);

% Set up fittype and options
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
ft = fittype(@(a,b,c,x) c + a*sin(2*pi*f*x + b),... % for a single fixed-frequency sinusoid
           'coefficients', {'a', 'b', 'c'});

% Estimate starting point values
fs = 1 / mean(diff(xData));
[~,fz,mag,phs] = chirpz(yData, xData, 0, fs/2);
[~,mI] = min(abs(fz - f));
mag_est = mag(mI);
phs_est = phs(mI);

% Bounds
opts.Lower = [-80 -pi -Inf];
opts.Upper = [80 pi Inf];
opts.StartPoint = [mag_est phs_est 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if showplot
    
    % Plot fit with data.
    %figure( 'Name', 'Fixed frequency sine fit' );
    figure (403)
    hold on
    title(gof.rsquare)
    h(1) = plot(x, y, 'k');
    h(2) = plot(xData, yData, 'b');
    h(3) = plot(fitresult, 'r');
    legend( h, 'raw', 'detrend', 'fit', 'Location', 'NorthEast', 'Interpreter', 'none' );

    % Label axes
    xlabel( 'time', 'Interpreter', 'none' );
    ylabel( 'y', 'Interpreter', 'none' );
    grid on
end

end


