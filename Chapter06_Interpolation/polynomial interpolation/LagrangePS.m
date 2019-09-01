% Lagrange Practical Solution to calculate the zero
% 10170437 Mark Taylor
function x0 = LagrangePS(table)

% swap y with x
table=table(:,[2 1]);

x0=LagrangeEvalAt(table,0);
end
