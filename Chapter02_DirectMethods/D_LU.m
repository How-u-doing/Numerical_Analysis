% Doolittle's LU decomposition: A=LU, where L is a unit lower triangular
% matrix, and U is an upper triangular matrix.
% 10170437 Mark Taylor
function [L,U]=D_LU(A)
% Implementation via Doolittle's Method without pivoting.

[m,n]=size(A);
if m ~= n
    error('A must be square!')
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
        else
            error('U(%d,%d)=0!\n',i,i)
        end
    end
end

% Combine L & U into one matrix if there is only one output argument.
if nargout==1
    L=tril(L,-1)+triu(U);
end

end


 