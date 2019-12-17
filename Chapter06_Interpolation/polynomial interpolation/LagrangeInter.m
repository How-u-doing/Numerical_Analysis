% Lagrange interpolation
function [s, L] = LagrangeInter(x,y,xq)
if nargin==3 && nargout<=1          % only ouput yq
    s = LagrangeNumer(x,y,xq);    	% numerically solve
    return
end
% input only with x & y or need to output the Lagrange polynomial
n=length(x);    % no. of sample points
syms t;
L=0;
for i=1:n
    Li=1;
    for k=1:n
        if k~=i
            Li=Li*(t-x(k))/(x(i)-x(k));
        end
    end
    L=L+Li*y(i);
end
L=simplify(L);
s=subs(L,t,xq);
s=double(s);


function s = LagrangeNumer(x, y, xq)
n=length(x);
s=0;
for i=1:n
    Li=1;
    for k=1:n
        if k~=i
            Li=Li.*(xq-x(k))/(x(i)-x(k));
        end
    end
    s=s+Li*y(i);
end

