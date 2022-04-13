function dx = sparseODE(t,x,A)
% Full structure expression of indetified ODES

d = size(A,1);   % number of rows, indicating the dim of ODEs

switch d
    case 1
        xLib = [ x(1)  x(1)*x(1)  x(1)*x(1)*x(1)]';
    case 2
        xLib = [ x(1)  x(1)*x(1)  x(1)*x(1)*x(1) ...  
                 x(2)  x(1)*x(2)  x(1)*x(1)*x(2)  x(2)*x(2)  x(1)*x(2)*x(2)  x(2)*x(2)*x(2)]';
    case 3
        xLib = [ x(1)  x(1)*x(1)  x(1)*x(1)*x(1) ...  
                 x(2)  x(1)*x(2)  x(1)*x(1)*x(2)  x(2)*x(2)  x(1)*x(2)*x(2)  x(2)*x(2)*x(2) ...
                 x(3)  x(1)*x(3)  x(1)*x(1)*x(3)  x(2)*x(3)  x(1)*x(2)*x(3)  x(2)*x(2)*x(3)  x(3)*x(3)  x(1)*x(3)*x(3)  x(2)*x(3)*x(3)  x(3)*x(3)*x(3) ]';
end


dx = A*xLib;

end

