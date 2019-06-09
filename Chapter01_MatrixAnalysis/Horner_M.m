% Matrix form of Horner's Method 
% 10170437 Mark Taylor
function Y = Horner_M(p,X) 

% Y = P(1)*X^N + P(2)*X^(N-1) + ... + P(N)*X + P(N+1)

[m,n] = size(X);
if m ~= n 
    error('Input must be square!')
end

np = length(p);
Y = zeros(m);
for i = 1:np
    Y = X * Y + diag(p(i) * ones(m,1));
end
