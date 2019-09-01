% ------------Newton Interpolation---------------
% Find a polynomial function of order "order" that goes through a set of points that are in the n-by-2 table
% 10170437 Mark Taylor
function fun = Newton(table,order)
m=size(table,1);            % the number of known points
if nargin<2
    order=m-1;
elseif order<1
    error('order must be no less than 1')
elseif order>=m
    error('order must be less than the number of given points')
end

n=2+order;
table(:,3:n)=zeros(m,order);

syms x;
N=sym('N',[order,1]);
fun=table(1,2);

j=3;
while j<=n
    i=j-1;
    while i<=m
    	table(i,j)=(table(i,j-1)-table(i-1,j-1))/(table(i,1)-table(i-(j-2),1));    
    	i=i+1;
    end
        
    N(j-2)=1;
    k=1;
    while k<j-1
       N(j-2)=N(j-2)*(x-table(k,1));
       k=k+1;
    end   
    fun=fun+table(j-1,j)*N(j-2);
    
    j=j+1;
end

%fun=matlabFunction(fun);
end

