% Row reduction with partial pivoting
% 10170437 Mark Taylor

function R = RowReduce(A)

[m,n]=size(A);
R=A;
i=1;
j=1;

while i<=m && j<=n
    MAX=abs(R(i,j));
    maxrn=i;                        % max row number    
    for k=i+1:m        
        if abs(R(k,j))>MAX
            MAX=abs(R(k,j));
            maxrn=k;
        end
    end
    
    if MAX<eps
       j=j+1;
       continue    
    elseif i~=maxrn
        temp=R(i,j:n);
        R(i,j:n)=R(maxrn,j:n);
        R(maxrn,j:n)=temp;
    end
        
    % Perform Gaussian elimination
    for t=i+1:m
        if abs(R(t,j))>eps
            R(t,j:n)=R(t,j:n)-R(t,j)/R(i,j)*R(i,j:n);
        end
    end
    i=i+1;
    j=j+1;
end

end