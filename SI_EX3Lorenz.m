%
clc
clear
close all

addpath('E:\fdaM') % "http://www.psych.mcgill.ca/misc/fda/" for basis functions
addpath('./odes')
addpath('./utils')

%% Sate-space equation

% STAGE I: state equation  
x0 = [-5; 10; 30];                	% initial vector 
sigma = 10; beta = 8/3; rho = 28;	% strauctral parameter (wiki)

nstat = length(x0);                 % dimension of state varibale

dt = 0.005; 
tobs = (0:dt:5)';                  % sampling time instant
nobs = length(tobs); 
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,nstat));
[~, xtru] = ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tobs,x0,options);

% STAGE II: observation equation 
for nvr = [5, 10, 15]/100           % noise level 5%, 10%, 15% 
    rng(1)                              % for reproducibility
    nois = nvr*std(xtru).*randn(size(xtru)); 
    xobs = xtru + nois;                 % noisy observations

    %% Identification
    pord = 3;
    lambda = 0.8;

    % noise-free case
    [Theta,~] = candLibrary(xtru,pord);
    xi0 = stlsIntg(xtru,Theta,lambda,dt);

    % noisy case
    [Theta,~] = candLibrary(xobs,pord);
    xi1 = stlsIntg(xobs,Theta,lambda,dt);

    % smoothing case
    for istat=1:nstat
        [xsmo(:,istat), ~] = pss_gcv(tobs,xobs(:,istat),0);
    end

    [Theta,pad] = candLibrary(xsmo,pord);
    xi2 = stlsIntg(xsmo,Theta,lambda,dt);

    array2table([ [nan(1,size(pad,2));pad], xi0 xi1 xi2], ...
                'VariableNames',{'x1_order','x2_order', 'x3_order',...
                                 'x1_coef_nf','x2_coef_nf','x3_coef_nf', ...
                                 'x1_coef_n','x2_coef_n', 'x3_coef_n', ...
                                 'x1_coef_sm','x2_coef_sm','x3_coef_sm'})
    %% Solution
    eta = xi2(1,:); 
    Xi = xi2(2:end,:)';
    tfor = [0:dt:15]';                  % whole time range
    [~, xfit] = ode45(@(t,x)sparseODE(t,x,Xi),tfor,eta,options);
    [~, xtrufor] = ode45(@(t,x)lorenz(t,x,sigma,beta,rho),tfor,x0,options);  

    %% Figure: phase
    f = figure;
    plot3(xobs(:,1),xobs(:,2),xobs(:,3),'g.','markersize',12); hold on
    plot3(xtrufor(:,1),xtrufor(:,2),xtrufor(:,3),'r-','linewidth',2.0)
    plot3(xfit(:,1),xfit(:,2),xfit(:,3),'k--','linewidth',2.0); hold off
    view(27,16); grid on; grid minor
    xlabel('$x_1(t)$','interpreter','latex')
    ylabel('$x_2(t)$','interpreter','latex')
    zlabel('$x_3(t)$','interpreter','latex')
    title([['nvr = ',num2str(nvr*100)],'%'])
    set(gca,'fontsize',13)
    set(gcf,'position',[100 200 500 500])

    set(f,'PaperSize',[15 10])
    print('-painters',f,['Lorenz', num2str(nvr*100)],'-dpdf')
end

%% Figure: time-series 
ff = figure;
for istat=1:nstat
    subplot(3,1,istat)    
    plot(tobs,xobs(:,istat),'g.','markersize',10); hold on
    plot(tfor,xtrufor(:,istat),'-r','linewidth',2.0); 
    plot(tfor,xfit(:,istat),'--k','linewidth',2.0); hold off
    switch istat
        case 1
            ylabel('$x_1(t)$','interpreter','latex')            
        case 2
            ylabel('$x_2(t)$','interpreter','latex')
        case 3 
            xlabel('$t$','interpreter','latex')
            ylabel('$x_3(t)$','interpreter','latex')
    end
    grid on; grid minor
    xline(5,'--b',{'Forecasting','horizon'},'LineWidth',2);
    set(gca,'fontsize',12)
    set(gcf,'position',[100 200 1000 500])
end
set(ff,'PaperSize',[15 10])
print(ff, ['Fore_Lorenz',num2str(nvr*100)],'-dpdf')






