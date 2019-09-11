% ------------Hermite Interpolation---------------
% Find a polynomial H(x) of at most£¨2m-1)th degree that goes through a set of points  
% that are in the m-by-3 table and with the same  derivatives as f(x) at those points.
% 10170437 Mark Taylor
function [fun, newTable] = Hermite(table)
% table=[X,f(X),f'(X)], where X, f(X) and f'(X) are given column vectors of m-th dimension

m=size(table,1);            % the number of known points
newTable=NaN(2*m,2*m+1);

% Initialization
for i=1:m
    newTable([2*i-1;2*i],[1,2])=table([i;i],[1,2]);
    newTable(2*i,3)=table(i,3);
end

[m,n]=size(newTable);
for i=3:2:m-1
    newTable(i,3)=(newTable(i,2)-newTable(i-1,2))/(newTable(i,1)-newTable(i-1,1));    
end
syms x;
fun=newTable(1,2)+newTable(2,3)*(x-newTable(1,1));

for j=4:n
    for i=j-1:m
    	newTable(i,j)=(newTable(i,j-1)-newTable(i-1,j-1))/(newTable(i,1)-newTable(i-(j-2),1));
    end
        
    N=1;
    for k=1:j-2
       N=N*(x-newTable(k,1));
    end
    
    fun=fun+newTable(j-1,j)*N;
end

%fun=simplify(fun);
fun=vpa(fun,7);
end
