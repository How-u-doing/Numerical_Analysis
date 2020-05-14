% Crout's LU Factorizaton: A=LU, where L is a lower triangular
% matrix and U is a unit upper triangular matrix
% 10170437 Mark Taylor
function [L,U]=C_LU(A) 
% Implementation via Crout's Method without pivoting.
    
[m,n] = size(A);
if m ~= n
    error('A must be square!')
end

U=eye(n);
L=zeros(n);

for i=1:n
    for k=i:n  
        L(k,i)=A(k,i)-L(k,1:i-1)*U(1:i-1,i);
    end
    for k=i+1:n
        if abs(L(i,i))>eps            
            U(i,k)=(A(i,k)-L(i,1:i-1)*U(1:i-1,k))/L(i,i);        
        else
            error('L(%d,%d)=0!\n',i,i)
        end
    end
end

% Combine L & U into one matrix if there is only one output argument.
if nargout==1
    L=tril(L)+triu(U,1);
end

end


 