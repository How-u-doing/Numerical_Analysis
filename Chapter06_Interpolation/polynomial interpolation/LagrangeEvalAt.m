% Evaluate the value at x of a Lagrange interpolating polynomial 
% which is determined by a set of points in arg table.
% 10170437 Mark Taylor
function key = LagrangeEvalAt(table,x)
% table=[X,f(X)], where X and f(X) are given column vectors of n-th dimension
% x: a scalar or a vector to be evaluated

n=size(table,1);            % the number of known points
key=0;
for i=1:n
    Li=1;
    for k=1:n
        if k==i
            continue            
        end
        % here use .* to combine the vector input x
        Li=Li.*(x-table(k,1))/(table(i,1)-table(k,1));
    end
        
    key=key+Li*table(i,2);
end

end

