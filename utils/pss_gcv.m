function [xs,lamb] = pss_gcv(t,yobs,gcv_plot)
% Written by Baolei Wei
%   This function outputs the spline smoothing values from the spline
%   interpolation, xs; the value of smoothing hyper-parameter, rho, is
%   determined by generalised cross-validation criteria. 

% Input: the time index vector (nx1) t;
%        the noisy time series vector (nx1) yobs;
% Output: the estimated time series vector (nx1) xs;
%         the optimal hyper-parameter, lamb = rho/(1-rho);

if nargin < 3
    gcv_plot = 0; 
end

l = -6; u = 4;
lamb = logspace(l,u,10*(u-l))';
gcv = zeros(length(lamb),1);
xs = zeros(length(t),length(lamb));

for ind=1:length(lamb)
    [xs(:,ind), gcv(ind)] = pps(t,yobs,lamb(ind));
end

[~,indmax] = min(gcv);
xs = xs(:,indmax);

if gcv_plot == 1
    figure('name','gcv')
    semilogx(lamb, gcv, '-b.','linewidth',1.0,'markersize',10)
%     title('GCV error','fontsize', 15)
    grid on; grid minor
    xlabel('\rho/(1-\rho)','fontsize',15)
    ylabel('GCV(\rho)','fontsize',15)
    set(gca,'fontname','book antiqua','fontsize',15)
    set(gcf,'position',[100 200 450 450])
end

end


function [xhat,gcv] = pps(t,yobs,rho)

%% basis functions
knots = t;
norder = 4;         % order of b-splines: 4 gives cubic splines
nobs = length(t);
nbasis = nobs+norder-2;     % number of basis functions
fbasis = create_bspline_basis([min(t),max(t)],nbasis,norder,knots);

%% calculate penalty term using Simpson's formula
% knots for quadrature. NOTE: the number of knots can be different from
% the numbers of time knots but need to be odd, here, 2*length(t)+1.
qnots = linspace(min(t),max(t),2*length(t)+1)'; 
% weights for Simpson's formula 
wnots = ones(length(qnots),1); 
wnots(2:2:end-1) = 4; 
wnots(3:2:end-2) = 2;
h = qnots(2)-qnots(1); 
wnots = wnots.*h/3;
% 2nd order derivatives of basis functions
D2Qmat = eval_basis(qnots, fbasis, 2);
Qmat   = D2Qmat'*(D2Qmat.*(wnots*ones(1,nbasis)));

%% estimates of state variable form observations
Rmat = eval_basis(t, fbasis);
Invmat = (Rmat'*Rmat+rho*Qmat)\eye(size(Rmat,2));
bhat = Invmat*(Rmat'*yobs);
xhat = Rmat*bhat;

%% generalized cross validation error
Smat = Rmat*Invmat*Rmat'; 
gcv1 = norm( (eye(nobs)-Smat)*yobs,2 )^2;
gcv2 = trace( eye(nobs)-Smat )^2;
gcv = nobs*gcv1/gcv2;

end





