%
clc
clear
close all

addpath('E:\fdaM') % "http://www.psych.mcgill.ca/misc/fda/" for basis functions
addpath('./odes')
addpath('./utils')

%% Sate-space equation

% STAGE I: state equation 
x0 = [1.8; 1.8];            % initial value & 
Be = [2/3 -4/3; -1 1];	% strauctral parameter (wiki)

nstat = length(x0);         % dimension of state varibale

h = 0.01; 
tobs = (0:h:10)';           % sampling time instant
nobs = length(tobs); 
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,nstat));
[~, xtru] = ode23(@(t,x)lotkav(t,x,Be),tobs,x0,options);

% STAGE II: observation equation 
for nvr = [5, 10, 15]/100          % .05,.10,.15
    rng(1)                      % for reproducibility
    nois = nvr*std(xtru).*randn(size(xtru)); 
    xobs = xtru + nois;         % noisy observations

    %% Structure Identification
    pord = 3;
    lambda = 0.3;

    % noise-free case
    [Theta,~] = candLibrary(xtru,pord);
    xi0 = stlsIntg(xtru,Theta,lambda,h);

    % noisy case
    [Theta,~] = candLibrary(xobs,pord);
    xi1 = stlsIntg(xobs,Theta,lambda,h);

    % smoothing case
    for istat=1:nstat
        [xsmo(:,istat), ~] = pss_gcv(tobs,xobs(:,istat),0);
    end

    [Theta,pad] = candLibrary(xsmo,pord);
    xi2 = stlsIntg(xsmo,Theta,lambda,h);

    array2table([ [nan(1,size(pad,2));pad], xi0 xi1 xi2], ...
                'VariableNames',{'x1_order','x2_order', 'x1_coef_nf','x2_coef_nf', ...
                                                        'x1_coef_n','x2_coef_n','x1_coef_sm','x2_coef_sm'})

    %% Trajectories
    % % PHASE
    x0 = xi2(1,:); 
    A = xi2(2:end,:)';
    [~, xfit] = ode23(@(t,x)sparseODE(t,x,A),tobs,x0,options);

    f = figure('name','trajectory');
    plot(xobs(:,1),xobs(:,2),'.g','markersize',10); hold on
    plot(xtru(:,1),xtru(:,2),'-r','linewidth',2.0)
    plot(xfit(:,1),xfit(:,2),'--k','linewidth',2.0); hold off
    grid on; grid minor
    xlabel('$x_1(t)$','interpreter','latex')
    ylabel('$x_2(t)$','interpreter','latex')
    title([['nvr = ',num2str(nvr*100)],'%'])
    set(gca,'fontsize',13)
    set(gcf,'position',[100 200 500 500])

    set(f,'PaperSize',[15 10])
    print(f, ['Lokta', num2str(nvr*100)],'-dpdf')
end

%% FIGURE
tfor = (0:h:30)';
[~, xtru] = ode23(@(t,x)lotkav(t,x,Be),tfor,x0,options);
[~, xfor] = ode23(@(t,x)sparseODE(t,x,A),tfor,x0,options);

ff = figure('name','states');
for istat=1:nstat
    subplot(2,1,istat)
    plot(tobs,xobs(:,istat),'g.','markersize',10); hold on
    plot(tfor,xtru(:,istat),'-r','linewidth',2.0); 
    plot(tfor,xfor(:,istat),'--k','linewidth',1.5); hold off 
    switch istat
        case 1
            ylabel('$x_1(t)$','interpreter','latex')
        case 2
            xlabel('$t$','interpreter','latex')
            ylabel('$x_2(t)$','interpreter','latex')
    end
    grid on; grid minor
    xline(10,'--b',{'Forecasting','horizon'},'LineWidth',2);
    set(gca,'fontsize',12)
    set(gcf,'position',[100 200 1000 500])
end
set(ff,'PaperSize',[15 10])
print(ff, ['Fore_Lokta',num2str(nvr*100)],'-dpdf')






