function [mdl] = taufxn(x,y)

%%% specify model function
model = @(a,x) a(1).*exp(-x(:,1)/a(2))+a(3); % exponential a*exp(-x/b)+c

% starting points for coeff. estimates
beta0 = [.01 1 -10]; 

mdl = fitnlm(x,y,model,beta0)

end