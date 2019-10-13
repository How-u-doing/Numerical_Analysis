% Bernstein Base functioon on interval [0,1]
% 10170437 Mark Taylor
function P=BernsteinBase(f,n,x)
P=f(0)*(1-x).^n;
for k=1:n
    P=P+f(k/n)*nchoosek(n,k)*x.^k.*(1-x).^(n-k);
end
end

