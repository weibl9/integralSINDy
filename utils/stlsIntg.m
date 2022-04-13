function xi = stlsIntg(x,Theta,lambda,dt)
% STLS for sparsity: 
%    lambda is the sparsification knob.

nstat = size(x,2);            % number of state variables 
[nobs,ncand] = size(Theta);   % number of candidate functions

% Cusum generation: Euler formula
Omega = nan(nobs-1,ncand);
for indcol = 1:ncand
    omega = Theta(:,indcol);
    Omega(:,indcol) = cumsum(omega(1:end-1)+omega(2:end))*dt/2;
end
Omega = [ones(nobs-1,1) Omega];

% sequential threshold least squares
xi = Omega\x(2:end,:);          % least-square initialization
for istat = 1:nstat             % updates of structural coefficients
    for k=1:10
        biginds = (abs(xi(:,istat))>lambda);
        biginds(1,:) = 1;                  % excluding initial value 
        xi(~biginds,istat) = 0;            % set small cofficients to 0
        
        xi(biginds,istat) = Omega(:,biginds)\x(2:end,istat); 
    end
end

end




