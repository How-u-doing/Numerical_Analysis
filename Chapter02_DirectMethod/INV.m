% Compute the inverse of square matrix A. 
% 10170437 Mark Taylor
function Inv_A = INV(A)

if nargin==0
    error('Please input a square matrix as the input to compute its inverse matrix!')    
end
    
[m,n]=size(A);
if m ~= n % make sure A is a square matrix to continue following steps
    error('Input error! Input must be a square matrix!')
end
U=[A,eye(n)];
Inv_A=zeros(n);
for j=1:m    
    % Search down from U(j:m,j) and find the entry in the left column with 
    % the largest absolute value as pivot,indicated as B(k,j)(k>=j), then 
    % swap row j with row k if k>j
    k=maxIndex(U(:,j),j,m);
    if U(k,j)~=0
        if k~=j
            temp=U(j,:);
            U(j,:)=U(k,:);
            U(k,:)=temp;               
        end 
    else 
        disp('Input is a singular matrix, which has no invese matrix!')
        return
    end   
    
%     Let first nonzero entry of each row be equal to one   
    U(j,:)=U(j,:)/U(j,j);
    
%     eliminate the entries that above U(j,j)
    if j>1 
        for t=1:j-1
            if U(t,j)~=0
                U(t,:)=U(t,:)-U(t,j)*U(j,:);
            end
        end
    end  
    
    
    for i=j+1:m 
        if(U(i,j)~=0)
            U(i,:)=U(i,:)-U(i,j)*U(j,:);
        end                 
    end
    
end


Inv_A=U(:,n+1:2*n);
end


