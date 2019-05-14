% Gaussian Elimination:  solve linear systems of algebraic 
% equations of the form Ax=b, where A is a nonsigular matrix.
% 10170437 Mark Taylor
function [U,x] = GauEli(A, b)
    
[m,n]=size(A);
if m ~= n % make sure A is a square matrix to continue following steps
    error('Input error! Input must be a square matrix!')
end
U=[A,b];% Augemented matrix of A

for j=1:n-1    
    % Search down from U(j:m,j) and find the entry in the left column with 
    % the largest absolute value as pivot,indicated as B(k,j)(k>=j), then 
    % swap row j with row k if k>j
    k=maxIndex(U(:,j),j,m);
    if U(k,j)~=0
        if k~=j
            temp=U(j,j:n+1);
            U(j,j:n+1)=U(k,j:n+1);
            U(k,j:n+1)=temp;               
        end 
    else 
        error('Input is a singular matrix!')
    end
    
    for i=j+1:m
        if(U(i,j)~=0)
            t=U(i,j)/U(j,j);
            U(i,j:n+1)=U(i,j:n+1)-t*U(j,j:n+1);
        end            
    end
end

x=zeros(m,1);
x(n)=U(m,n+1)/U(m,n);
for i=n-1:-1:1
    x(i)=(U(i,n+1)-U(i,i+1:n)*x(i+1:n))/U(i,i);
end

end


