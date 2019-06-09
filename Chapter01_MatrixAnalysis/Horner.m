% Horner's Method
% 10170437 Mark Taylor
function y = Horner(p,x)  

% y = P(1)*X^N + P(2)*X^(N-1) + ... + P(N)*X + P(N+1)


np=length(p);
y=zeros(size(x));
for k=1:np
    y=x.*y+p(k);
end
end
