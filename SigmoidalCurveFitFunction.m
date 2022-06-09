function y=SigmoidalCurveFitFunction(x,a,b,c,n)
y= a.*(((x+20).^n)./(((x+20).^n)+(b.^n)))+c;
end