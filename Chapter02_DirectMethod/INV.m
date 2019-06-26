% Compute the inverse of square matrix A. 
% 10170437 Mark Taylor
function Inv_A = INV(A)

[m,n]=size(A);
if m ~= n
    error('A must be square!')
end
U=[A,eye(n)];
n1=2*n;
for j=1:n
    k=maxIndex(U(:,j),j,n);
    if abs(U(k,j))>eps
        if k~=j
            temp=U(j,j:n1);
            U(j,j:n1)=U(k,j:n1);
            U(k,j:n1)=temp;               
        end 
    else 
        Inv_A='A is singular!';
        return
    end   
    
    % Render first nonzero entry of each row to be one.
    U(j,:)=U(j,:)/U(j,j);
    
    % Eliminate the entries that above U(j,j).
    if j>1
        for t=1:j-1
            if abs(U(t,j))>eps
                U(t,j:n1)=U(t,j:n1)-U(t,j)*U(j,j:n1);
            end
        end
    end  
    
    % Eliminate the entries that below U(j,j).
    for i=j+1:n 
        if abs(U(i,j))>eps
            U(i,j:n1)=U(i,j:n1)-U(i,j)*U(j,j:n1);
        end                 
    end
    
end


Inv_A=U(:,n+1:n1);
end