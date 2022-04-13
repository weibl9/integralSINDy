function [Omega,pad] = candLibrary(X,pord) 
	% Generate candidate librady. 
	% X:  multi-variate time series of state variable
	% X = [x1 x2 x3 ... xnvar] 
	% pord: ploynomial order for basis expansion (integer)

    [~,nvar] = size(X);    
    
	%% polynomial basis expansion
    if pord == 1
        pad = [1:pord]';            % index combination
        Omega = x2fx(X,pad);        % generate design matrix
    else
        [iniPerm{1:nvar}] = ndgrid(0:pord);
        allPerm = reshape(cat(nvar+1,iniPerm{:}), [],nvar);   % permutations nvar^(pr+1)

        ind = sum(allPerm,2)<=pord;                           % satisfy condition
        ind(1) = 0;                                           % 1<= power(pord) <= pr 

        pad = allPerm(ind,:);                                 % all candidate perms
        Omega = x2fx(X,pad);                                  % generate design matrix
    end

end




















