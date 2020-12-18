% Gaussian Elimination: solve linear systems 
% of algebraic equations of the form Ax=b.
% 10170437 Mark Taylor

function [x,U] = GauEli(A, b)

[m,n]=size(A);
if m ~= n 
    error('A must be a square matrix!')
end

U=[A,b];        % Augemented matrix of A

for j=1:n
    k=maxIndex(U(:,j),j,n);
    if abs(U(k,j))>eps
        if k~=j
            temp=U(j,j:n+1);
            U(j,j:n+1)=U(k,j:n+1);
            U(k,j:n+1)=temp;               
        end 
    else 
        x='The system has infinitely many solutions!';
        U=NaN;
        return;        
    end
    
    % Perform Gauss elimination.
    for i=j+1:n
        if abs(U(i,j))>eps
            t=U(i,j)/U(j,j);
            U(i,j:n+1)=U(i,j:n+1)-t*U(j,j:n+1);
        end
    end
end

x=zeros(n,1);
% Sovle upper diagonal linear system, i.e. U(1:n,1:n)x=U(:,n+1).
x(n)=U(n,n+1)/U(n,n);
for i=n-1:-1:1
    x(i)=(U(i,n+1)-U(i,i+1:n)*x(i+1:n))/U(i,i);
end

end