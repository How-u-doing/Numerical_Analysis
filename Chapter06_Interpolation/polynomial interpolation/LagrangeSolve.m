% Give a set of known points of a function to find the approximate zero(s).
% Implementation by constructing a Lagrange interpolating polynomial
%ги x=P(y), calculate the function value at y=0 )
% See also LagrangePS.m which substitutes syms x with the exact number x
% in function LagrangeEvalAt(), and which consequently computes faster.
% 10170437 Mark Taylor
function x0 = LagrangeSolve(table)

% swap y with x
table=table(:,[2 1]);
L=Lagrange(table);

% following  <->  x0=subs(L,0), which usually produces a symbolic fraction,
% e.g. 1247/1250. Using x0=double(subs(L,0)) gives rise to a deciamal  
L=matlabFunction(L);
x0=L(0); % compute the value at y=0, generating a decimal

end
