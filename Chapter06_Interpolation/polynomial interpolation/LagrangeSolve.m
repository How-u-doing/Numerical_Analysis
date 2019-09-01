% Give a set of known points of a function to find the approximate zero(s).
% Implementation by constructing a Lagrange interpolating polynomial
%£¨x=P(y), calculate the function value at y=0)
% See also LagrangePS.m for a more practical & effective solution.
% 10170437 Mark Taylor
function x0 = LagrangeSolve(table)

% swap y with x
table=table(:,[2 1]);

L=Lagrange(table);
L=matlabFunction(L);
x0=L(0); % compute the value at y=0
end
