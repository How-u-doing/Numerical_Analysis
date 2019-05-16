% Solve Ux=b, where U is an upper triangular matrix
% 10170437 Mark Taylor
function x=solveTriu(U,b)
    [m,n]=size(U);
    if m~=n % For implementation efficiency, we do not check if U is upper triangular
        error('First input must be a lower triangular matirx!')
    elseif length(b)~=n
        error('Dimension mismatch between column vector b and U!')
        
    end
    
    if all(diag(U))==false
        error('U is singular, and resultant no solution or infinite solutions in Ux=b!')
    end 
    
    x=zeros(n,1);
    x(n)=b(n)/U(n,n);
    for i=n-1:-1:1
        x(i)=(b(i)-U(i,i+1:n)*x(i+1:n))/U(i,i);
    end
    
end    