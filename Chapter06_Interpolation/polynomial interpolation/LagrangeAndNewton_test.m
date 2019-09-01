% Lagrange and Newton test
% 10170437 Mark Taylor

fprintf('function Lagrange() & Newton() test:\n');
table=[-3,-1;0,2;3,2;6,10]
L=Lagrange(table)
N=Newton(table)


fprintf('\n\nfunction LagrangeSolve() & LagrangePS() test:\n');
% f(x)=x^4-3x+1
table2=[0.1,0.70010;0.2,0.40160;0.3,0.10810;0.4,-0.17440]
% If you run these two functions in command window respectively£¨
% you'll find the latter runs much faster than the former one.
% Reason's simple: latter doesn't need to call syms family, which
% saves a lot of time.
x0=LagrangeSolve(table2)
x1=LagrangePS(table2)


