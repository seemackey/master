function [MyCoeffs] = taufxn_v2(x,y,tauguess,rangeguess)


% inputs are the data and a guess for 

xData=x';
yData=y;

% settings for the fit
ft = fittype( 'a*exp(-x/b)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [rangeguess tauguess -18];
opts.Upper = [100 150 60];
opts.StartPoint = [0.01 1 -1];



% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
clf;
A = fitresult;

% gets coefficients to calculate dynamic range of psych fxn
MyCoeffs = coeffvalues(fitresult);
end