% ------------Newton Interpolation---------------
% Find a polynomial of order-th order that goes through a set of points that are in the m-by-2 table
% 10170437 Mark Taylor
function [fun,table] = Newton(table,order)
% table=[X,f(X)], where X and f(X) are given column vectors of m-th dimension
% order is the order of the Newton interpolating polynomial to be returned

m=size(table,1);            % the number of known points
if nargin<2
    order=m-1;              % to calculate the Newton interpolating polynomial of (m-1)-th order by default
elseif order<1
    % when order<m-1, use preceding (order+1) points to approximate the polynomial
    error('order must be no less than 1')
elseif order>=m
    error('order must be less than the number of given points')
end

n=2+order;
table(:,3:n)=NaN(m,order);

syms x;
fun=table(1,2);
for j=3:n
    for i=j-1:order+1
    	table(i,j)=(table(i,j-1)-table(i-1,j-1))/(table(i,1)-table(i-(j-2),1));
    end
        
    N=1;
    for k=1:j-2
       N=N*(x-table(k,1));
    end
    
    fun=fun+table(j-1,j)*N;
end

%fun=matlabFunction(fun);
end

