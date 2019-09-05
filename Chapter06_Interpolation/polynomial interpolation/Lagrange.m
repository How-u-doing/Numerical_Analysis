% ------------Lagrange Interpolation---------------
% Find a polynomial that goes through a set of points that are in the n-by-2 table
% 10170437 Mark Taylor
function fun = Lagrange(table)
% table=[X,f(X)], where X and f(X) are given column vectors of n-th dimension

n=size(table,1);            % the number of known points
syms x;
fun=0;
for i=1:n
    denominator=1; 
    Li=1;
    for k=1:n        
        if k==i
            continue            
        end
        
        % Here we split numerator  and  denominator  in  order  to make the 
        % coefficient of each item concerning x be one ( i.e. x-xi , rather 
        % than  x/3 - xi/3  or so, you will see it in the final expression.
        % This has something to do with the string handling mechanics built
        % in MATLAB )
        Li=Li*(x-table(k,1));
        denominator=denominator*(table(i,1)-table(k,1));
    end
    
    fun=fun+Li/denominator*table(i,2);
end

%fun=simplify(fun);
%fun=matlabFunction(fun);
end

