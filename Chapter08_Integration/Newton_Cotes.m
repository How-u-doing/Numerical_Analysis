% Newton-Cotes integration formula
function s = Newton_Cotes(f,a,b,n)
h=(b-a)/n;
x=a:h:b;
% construct a Lagrange polynomial based on n+1 evenly spaced points(i.e. x)
syms t;
L=0;
for i=1:n+1
    Li=1;
    for k=1:n+1
        if k~=i
            Li=Li*(t-x(k))/(x(i)-x(k));
        end
    end
    L=L+Li*f(x(i));
end
% integrate Lagrange polynomial
s=int(L,a,b);
s=double(s);    % return numerical integral
