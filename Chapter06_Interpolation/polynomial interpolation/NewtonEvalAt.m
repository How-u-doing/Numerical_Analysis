% Evaluate the value at x of a Newton interpolating polynomial of degree-th
% degree which is determined by the preceding (degree+1) points in arg table.
% 10170437 Mark Taylor
function key = NewtonEvalAt(table,x,degree)
% table=[X,f(X)], where X and f(X) are given column vectors of m-th dimension
% x is the point to be evaluated
% degree is the degree of the Newton interpolating polynomial to be approximated as

m=size(table,1);            % the number of known points
if nargin<3
    degree=m-1;             % approximation via the Newton interpolating polynomial of (m-1)-th degree by default
elseif degree<1
    % when degree<m-1, use preceding (degree+1) points to approximate the polynomial
    error('degree must be no less than 1')
elseif degree>=m
    error('degree must be less than the number of given points')
end

n=2+degree;
table(:,3:n)=NaN(m,degree);

key=table(1,2);
for j=3:n
    for i=j-1:degree+1
    	table(i,j)=(table(i,j-1)-table(i-1,j-1))/(table(i,1)-table(i-(j-2),1));
    end
        
    N=1;
    for k=1:j-2
       N=N*(x-table(k,1));
    end
    
    key=key+table(j-1,j)*N;
end

end

