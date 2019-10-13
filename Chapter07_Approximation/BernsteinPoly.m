% Bernstein Polynomial on interval [a,b]
% 10170437 Mark Taylor
function P=BernsteinPoly(f,n,x,a,b)
if nargin==3
    a=0;b=1;
end
% take a linear transform: [a,b]->[0,1]£¬f(x)->g(t) (= f(a+(b-a)t))
t=(x-a)/(b-a);     % t belongs to [0,1]
P=f(a)*(1-t).^n;
for k=1:n
    P=P+f(a+k/n*(b-a))*nchoosek(n,k)*t.^k.*(1-t).^(n-k);
end
end

