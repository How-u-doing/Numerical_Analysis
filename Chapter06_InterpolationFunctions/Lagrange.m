% Find a polynomial function that goes through a set of points that is in n-by-2 table
% 10170437 Mark Taylor
function fun = Lagrange(table)
n=size(table,1);

syms x;
L=sym('L',[n,1]);
fun=0;
i=1;
while i<=n
    denominator=1; 
    L(i)=1;
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
        
    fun=fun+L(i)/denominator*table(i,2);
    i=i+1;
end

%fun=matlabFunction(fun);
end

