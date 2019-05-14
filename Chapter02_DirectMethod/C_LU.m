% Crout's LU decomposition: A=LU, where L is a lower triangular
% matrix and U is a unit upper triangular matrix
% 10170437 Mark Taylor
function [L,U]=C_LU(A) 
% Implementation via Crout's Method without selecting pivot
    
[m,n] = size(A);
if m ~= n % make sure A is a square matrix to continue following steps
    error('Input error! Input must be a square matrix!')
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
        else  % Cannot be implemented without selecting pivot
            Error=sprintf('\nError: L(%d,%d)=0!\nCrout''s Method cannot be implemented without selecting pivot',i,i);
            disp(Error);
            L=NaN;
            U=NaN;
            return;
        end
    end
end

% Combine L & U into one matrix if there is only one output argument
if nargout==1
    L=tril(L)+triu(U)-eye(n);% L=tril(L)+triu(U,1)
end

end


 
