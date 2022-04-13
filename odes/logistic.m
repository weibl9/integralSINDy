function dx = logistic(t,x,beta)
% Logistic ODE: dx = a*x + b*x^2
    dx = beta(1)*x + beta(2)*x.^2;
end
