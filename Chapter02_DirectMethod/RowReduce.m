% Row reduction method,selecting entries with largest absolute value as column pivot
% 10170437 Mark Taylor
function [U,x,rank,fvi] = RowReduce(A, b)
% U is an upper triangular matrix, rank=rank([A,b])
if nargin==0
    error('Error: input cannot be NULL!')   
end

[m,n]=size(A);
if nargin==2
    if size(b,2)==1 % check if b is a column vector 
        U=[A,b];% Augemented matrix of A 
    else 
        error('The seceond argument input must be a column vector,whose length equals to the first demension of first input argument')   
    end
elseif nargin==1 
    b=zeros(m,1);
    U=[A,b];
    n=n-1;
end

i=1;j=1;
rank=0;
freeVarIndex=0;
% Perform Gaussian elimination on the augemented matrix until it is in reduced row echelon form
while i<=m && j<=n   
    % Search down from U(j:m,j) and find the entry in the left column with 
    % the largest absolute value as pivot,indicated as B(k,j) (k>=j), then 
    % swap row j with row k if k>j.
    MAX=abs(U(i,j));
    max_rn=i; % max row number    
    for k=i+1:m        
        if abs(U(k,j))>MAX
            MAX=abs(U(k,j));
            max_rn=k;
        end
    end
    
    if MAX==0 % all entries below B(i,j)(include B(i,j)) is equal to zero
       freeVarIndex(end+1)=j; 
       j=j+1; 
       continue    
    elseif i~=max_rn % switch current row i with row max_rn
        temp=U(i,j:n+1);
        U(i,j:n+1)=U(max_rn,j:n+1);
        U(max_rn,j:n+1)=temp;
    end
    % let first nonzero entry of each row be equal to one,of course here row#<=rank(B)
    U(i,j:n+1)=U(i,j:n+1)/U(i,j);
    % Gaussian elimination
    for t=i+1:m
        if(U(t,j)~=0)
            U(t,j:n+1)=U(t,j:n+1)-U(t,j)*U(i,j:n+1);
        end            
    end
    i=i+1;
    j=j+1;
    rank=rank+1;
end

if nargin==2    
    fvi=freeVarIndex(2:end); 
    if m==n && rank==m % Solve this linear system Ax=b if A is nonsingular
        x=zeros(m,1);
        x(n)=U(m,n+1)/U(m,n);
        for i=n-1:-1:1
            x(i)=(U(i,n+1)-U(i,i+1:n)*x(i+1:n))/U(i,i);
        end
    else % exists infinite solutions
        x='This linear system Ax=b exists infinite solutions';
        fprintf('\n') 
    end
end

if nargin==1 % [R,rank,fvi]=RowReduce(A)
   U(:,size(U,2))=[];
   fvi=freeVarIndex(2:end-1);
   x=rank;
   rank=fvi;         
end

end

