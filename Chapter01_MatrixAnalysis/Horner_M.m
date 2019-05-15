% Horner Method for matrix
% 10170437 Mark Taylor
function Y = Qin_M(p,X) 

% Y = P(1)*X^N + P(2)*X^(N-1) + ... + P(N)*X + P(N+1)

[m,n] = size(X);
if m ~= n % make sure input is a square matrix to continue following steps
    error('Input must be a square matrix!')
end

np = length(p);
Y = zeros(m);
for i = 1:np
    Y = X * Y + diag(p(i) * ones(m,1));
end
