% Solve Ly=b, where L is lower triangular 
% 10170437 Mark Taylor
function y=solveTril(L,b)
    [m,n]=size(L);
    if m~=n % For implementation efficiency, we do not check if L is lower triangular
        error('First input must be a lower triangular matirx!')
    elseif length(b)~=n
        error('Dimension mismatch between column vector b and L!')
        
    end
    
    if all(diag(L))==false
        error('L is singular, and resultant no solution or infinite solutions in Ly=b!')
    end 
    
    y=zeros(n,1);
    y(1)=b(1)/L(1,1);
    for i=2:n
        y(i)=(b(i)-L(i,1:i-1)*y(1:i-1))/L(i,i);
    end
    
end    