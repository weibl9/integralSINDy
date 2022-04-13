function dx = lotkav(t, x, beta)
% Lotka-Volterra predator-prey model
% beta = [beta1 beta2; beta3 beta4];
% dx = diag([1 - .01*x(2), -1 + .02*x(1)])*x;
dx = [
        beta(1,1)*x(1)+beta(1,2)*x(1)*x(2);
        beta(2,1)*x(2)+beta(2,2)*x(1)*x(2);
     ];
end