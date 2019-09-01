% Evaluate the value at x of a Newton  interpolation  polynomial of 
% order "order" which is determined by a set of points in arg table.
% 10170437 Mark Taylor
function key = NewtonEvalAt(table,x,order)
m=size(table,1);            % the number of known points
if nargin<3
    order=m-1;
elseif order<1
    error('order must be no less than 1')
elseif order>=m
    error('order must be less than the number of given points')
end

n=2+order;
table(:,3:n)=zeros(m,order);


N=ones(order,1);
key=table(1,2);

j=3;
while j<=n
    i=j-1;
    while i<=m
    	table(i,j)=(table(i,j-1)-table(i-1,j-1))/(table(i,1)-table(i-(j-2),1));    
    	i=i+1;
    end
      
    k=1;
    while k<j-1
       N(j-2)=N(j-2)*(x-table(k,1));
       k=k+1;
    end   
    key=key+table(j-1,j)*N(j-2);
    
    j=j+1;
end

end

