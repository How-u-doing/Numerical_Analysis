% Hermite Interpolation--according to our text's algorithm
% 10170437 Mark Taylor
function [y0,fun] = HI(X,Y,Yd,x0)
% [X,Y,Yd]=[X,f(X),f'(X)], where X is a nth dimention vector
% x0: a scalar or a vector to be approximated at

n=length(X);
syms x;
fun=0;
for i=1:n
    s=0;
    t=1;
    for j=1:n
               
        if j~=i
            t=t*((x-X(j))/(X(i)-X(j)))^2;
            s=s+1/(X(i)-X(j));                    
        end
        
    end
    
    fun=fun+Y(i)*t*(1-2*s*(x-X(i)))+Yd(i)*t*(x-X(i));
end

fun=simplify(fun);
fun=collect(fun);
fun=vpa(fun,7);

if nargin==4
    y0=subs(fun,x,x0);
    y0=vpa(y0,7);
else % nargin==3
    y0=fun;
end

end
