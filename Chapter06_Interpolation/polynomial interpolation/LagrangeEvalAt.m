% Evaluate the value at x of a Lagrange interpolation polynomial which
% is determined by a set of points in arg table.
% 10170437 Mark Taylor
function key = LagrangeEvalAt(table,x)
n=size(table,1);            % the number of known points

L=ones(n,1);
key=0;
i=1;
while i<=n
    denominator=1; 
    k=1;
    while k<=n
        if k==i
            k=k+1;
            continue            
        end
        
        L(i)=L(i)*(x-table(k,1));
        denominator=denominator*(table(i,1)-table(k,1));
        
        k=k+1;
    end
        
    key=key+L(i)/denominator*table(i,2);
    i=i+1;
end

end

