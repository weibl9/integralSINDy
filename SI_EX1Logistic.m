%
clc 
clear 
close all

addpath('E:\fdaM') % "http://www.psych.mcgill.ca/misc/fda/" for basis functions
addpath('./odes')
addpath('./utils')

%% Sate-space equation 
format short

% STAGE I: state equation 
x0 = 0.1;                  % initial value
beta = [1.6 -1.0];         % strauctral parameter

h = 0.01;
tobs = (0:h:4.5)';           % sampling time instant 
nobs = length(tobs); 
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,1));
[~, xtru] = ode45(@(t,x)logistic(t,x,beta),tobs,x0,options);

% STAGE II: observation equation
for nvr = [.05,.1,.2,.3,.4,.5]
    % nvr = 0.2;                  % 0.1--0.5
    rng(1)                      % for reproducibility
    nois = nvr*std(xtru)*randn(size(xtru)); 
    xobs = xtru + nois;         % noisy observations

    %% Structure Identification
    pord = 3;
    lambda = 0.1;

    % noise-free case
    [Theta,~] = candLibrary(xtru,pord);
    xi0 = stlsIntg(xtru,Theta,lambda,h);
%     array2table([ [nan(1,size(pad,2));pad] xi0], 'VariableNames',{'x_order','x_coef'})

    % noisy case
    [Theta,~] = candLibrary(xobs,pord);
    xi1 = stlsIntg(xobs,Theta,lambda,h);
%     array2table([ [nan(1,size(pad,2));pad] xi1], 'VariableNames',{'x_order','x_coef'})

    % smoothing case
    [xsmo, ~] = pss_gcv(tobs,xobs,0);
    [Theta,pad] = candLibrary(xsmo,pord);
    xi2 = stlsIntg(xsmo,Theta,lambda,h);
    
    array2table([ [nan(1,size(pad,2));pad] xi0 xi1 xi2], ...
        'VariableNames',{'x_order','x_coef_nf','x_coef_n','x_coef_sm'})

    %% Trajectories
    x0 = xi2(1); A = xi2(2:end)';
    [~, xfit] = ode45(@(t,x)sparseODE(t,x,A),tobs,x0,options);

    f = figure;
    plot(tobs,xobs,'.g','markersize',10); hold on
    plot(tobs,xtru,'-r','linewidth',2.0);
    plot(tobs,xfit,'--k','linewidth',2.0); hold off
    grid on; grid minor
    xlabel('$t$','interpreter','latex')
    ylabel('$x(t)$','interpreter','latex')
%     legend({'observations','true trajectory','identified trajectory'},'location','southeast')
    title([['nvr = ',num2str(nvr*100)],'%'])
    set(gca,'fontsize',13)
    set(gcf,'position',[100 200 450 450])
    
    set(f,'PaperSize',[15 10])
    print(f, ['Logistic', num2str(nvr*100)],'-dpdf')
    
    clear f
end

%% Figure
tfor = (0:h:6)';
[~, xtru] = ode45(@(t,x)sparseODE(t,x,A),tfor,x0,options);
[~, xfor] = ode45(@(t,x)sparseODE(t,x,A),tfor,x0,options);

ff = figure('name','states');
plot(tobs,xobs,'g.','markersize',10); hold on
plot(tfor,xtru,'-r','linewidth',2.0); 
plot(tfor,xfor,'--k','linewidth',1.5); hold off 
xlabel('$t$','interpreter','latex')
ylabel('$x(t)$','interpreter','latex')
grid on; grid minor
xline(4.5,'--b',{'Forecasting','horizon'},'LineWidth',2);
set(gca,'fontsize',12)
set(gcf,'position',[100 200 450 450])
set(ff,'PaperSize',[15 10])
print(ff, ['Fore_Logistic',num2str(nvr*100)],'-dpdf')









