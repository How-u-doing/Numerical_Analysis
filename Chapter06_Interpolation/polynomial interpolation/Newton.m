% ------------Newton Interpolation---------------
% Find a polynomial of degree-th degree that goes through a set of points that are in the m-by-2 table
% 10170437 Mark Taylor
function [fun,table] = Newton(table,degree)
% table=[X,f(X)], where X and f(X) are given column vectors of m-th dimension
% degree is the degree of the Newton interpolating polynomial to be returned

m=size(table,1);            % the number of known points
if nargin<2
    degree=m-1;              % to calculate the Newton interpolating polynomial of (m-1)-th degree by default
elseif degree<1
    % when degree<m-1, use preceding (degree+1) points to approximate the polynomial
    error('degree must be no less than 1')
elseif degree>=m
    error('degree must be less than the number of given points')
end

n=2+degree;
table(:,3:n)=NaN(m,degree);

syms x;
fun=table(1,2);
for j=3:n
    for i=j-1:degree+1
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

