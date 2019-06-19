% Gaussian Elimination: solve linear systems 
% of algebraic equations of the form Ax=b.
% 10170437 Mark Taylor

function [x,U] = GauEli(A, b)
% x is the unique solution(n-by-1 vector) or a message 
% that says the system has infinitely many solutions.

[m,n]=size(A);
if m ~= n 
    error('A must be a square matrix!')
end

U=[A,b];% Augemented matrix of A

for j=1:n-1 
    % Partial pivoting
    % Find the entry in the left column, i.e. U(j:m,j), with the largest 
    % absolute value, denoted by U(k,j)(k>=j), and set it as principal  
    % element, then swap row j with row k if k>j.
    k=maxIndex(U(:,j),j,m);
    if abs(U(k,j))>eps
        if k~=j
            temp=U(j,j:n+1);
            U(j,j:n+1)=U(k,j:n+1);
            U(k,j:n+1)=temp;               
        end 
    else 
        x='The system has infinitely many solutions!';
        return;        
    end
    
    % Perform Gauss elimination.
    for i=j+1:m
        if(U(i,j)~=0)
            t=U(i,j)/U(j,j);
            U(i,j:n+1)=U(i,j:n+1)-t*U(j,j:n+1);
        end            
    end
end

x=zeros(m,1);
% Sovle upper diagonal linear system, i.e. U(1:n,1:n)x=U(:,n+1).
x(n)=U(m,n+1)/U(m,n);
for i=n-1:-1:1
    x(i)=(U(i,n+1)-U(i,i+1:n)*x(i+1:n))/U(i,i);
end

end


