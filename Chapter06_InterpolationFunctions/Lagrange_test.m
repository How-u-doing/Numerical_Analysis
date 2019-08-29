% Lagrange_test & LagrangeSolve_test
% 10170437 Mark Taylor

fprintf('function Lagrange() test:\n');
table=[-3,-1;0,2;3,2;6,10]
L=Lagrange(table)


fprintf('\n\nfunction LagrangeSolve() test:\n');
% f(x)=x^4-3x+1
table2=[0.1,0.70010;0.2,0.40160;0.3,0.10810;0.4,-0.17440]
x0=LagrangeSolve(table2)


