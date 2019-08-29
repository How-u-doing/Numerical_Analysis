% Find the approximate solution by constructing a interpolation function
% 10170437 Mark Taylor
function x0 = LagrangeSolve(table)

table=table(:,[2 1]);
L=Lagrange(table);
L=matlabFunction(L);
x0=L(0); % compute the value at y=0
end
