% Doolittle's LU decomposition: A=LU, where L is a unit lower triangular
% matrix, and U is an upper triangular matrix.
% 10170437 Mark Taylor
function [L,U]=D_LU(A) 
% Implementation via Doolittle's Method without selecting pivot
    
[m,n] = size(A);
if m ~= n % make sure A is a square matrix to continue following steps
    error('Input error! Input must be a square matrix!')
end

L=eye(n);
U=zeros(n);

for i=1:n
    for k=i:n
        U(i,k)=A(i,k)-L(i,1:i-1)*U(1:i-1,k);        
    end
    for k=i+1:n
        if abs(U(i,i))>eps
            L(k,i)=(A(k,i)-L(k,1:i-1)*U(1:i-1,i))/U(i,i);
        else  % Cannot be implemented without selecting pivot
            Error=sprintf('\nError: U(%d,%d)=0!\nDoolittle''s Method cannot be implemented without selecting pivot',i,i);
            disp(Error);
            L=NaN;
            U=NaN;
            return;
        end
    end
end

% Combine L & U into one matrix as output if there is only one output argument
if nargout==1
    L=tril(L)-eye(n)+triu(U);% L=tril(L,-1)+triu(U)
end

end


 
