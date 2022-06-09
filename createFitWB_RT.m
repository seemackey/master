function [] = createFitWB_RT(HR,PC,threshold)
%CREATEFIT(LL,PC)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : LL
%      Y Output: PC
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 07-Jan-2019 23:14:35


%% Fit: 'untitled fit 1'.
% load('Array.mat')
% close figure 3
LLfit = HR(1:end,1)';
% minLev = min(LLfit);
% maxLev = max(LLfit);
% if minLev<0
%     LLfit = LLfit - minLev;
% end

PCfit = PC';
[xData, yData] = prepareCurveData( LLfit, PCfit );
a = threshold; % keeps threshold from being forced to fit
% Set up fittype and options.
ft = fittype( 'c-d.*exp(-(x./a)^b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [a 1 0.6 0.25];
opts.Upper = [600 50 1 1];
opts.StartPoint = [0 0.01 0.505508072942099 0.757740130578333];

figure( 'Name', 'Weibull CDF' );
Checkfit = 2;
% while Checkfit == 2


% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
clf;
A = fitresult;

% gets coefficients to calculate dynamic range of cdf
MyCoeffs = coeffvalues(fitresult);
Max=MyCoeffs(3)-MyCoeffs(4)*0.1; % max = saturation point - 10% of range, d
Min=MyCoeffs(3)-MyCoeffs(4)*0.9; % min = saturation point - 90% of range, d (coeff 4)

% computes dyn range of cdf
syms x;
eqn1=MyCoeffs(3)-MyCoeffs(4).*exp(-(x./MyCoeffs(1))^MyCoeffs(2))==(Max);
eqn2=MyCoeffs(3)-MyCoeffs(4).*exp(-(x./MyCoeffs(1))^MyCoeffs(2))==(Min);
 Xat90 = (solve(eqn1,x));
 Xat10 = (solve(eqn2,x));
DRat90 = abs(double(Xat90(1)));
DRat10 = abs(double(Xat10(1)));
DR=DRat90-DRat10;
RTcdfSlope=(Max-Min)/(DRat90-DRat10) % use dyn range to calc. slope


% if minLev<0
% Xax =  minLev:0.1:maxLev;
% 
% plot(Xax,eq)
% else

h = plot( fitresult, xData, yData );
legend( h, 'cdf vs. RT', 'Fit', 'Location', 'NorthEast' );
xlabel RT
ylabel cdf
grid on;

end
% Checkfit =input('Enter 1 if fit is good or Enter 2: ');

% end
