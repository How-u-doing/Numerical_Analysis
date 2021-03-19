% Exercise 2.7.2 in our text, p82
% $-\Delta{u}-k^2\,u=1$, in $(0,1)\times(0,1)$,
% $u=0$, on $\Gamma_1=\{x=0,0 \leq y \leq 1\}\cup\{0 \leq x \leq 1,y=1\}$
% $\nabla{u}\cdot\boldsymbol{n}=0$, on $\Gamma_2=\{0 \leq x \leq 1,y=0\}\cup\{x=1,0 \leq y \leq 1\}$

for k = [1,5,10,15,20,25]
    run_exercise_2_7_2(k)
end

%=================================================================

function run_exercise_2_7_2(k)
x=[0 1]; y=[0 1];
n_x=50; n_y=50;

[p,t] = generateMesh(x,y,n_x,n_y);

A = assembleMatrix(p,t,-k^2);
f = @(x,y) 1+x.*0+y.*0;
phi = assembleVector(p,t,f);

% process Dirichlet boundary condition
db_func = @(x,y) x.*0+y.*0;
[A,phi] = processDirichletBoundary(A,phi,p,db_func);

% solve the linear systems of equations
u = A\phi;

t=[t';zeros(1,size(t,1))];
p=p'; % required to be 2-by-nP

figure
pdemesh(p,0123,t,u) % 2nd arg `e` is unused, but has to be a matrix
title(['Helmholtz equation with k=',num2str(k)])
xlabel('x')
ylabel('y')
fig2svg(['../svg/Helmholtz_k=',num2str(k),'.svg']);

end

function [A,phi] = processDirichletBoundary(A,phi,p,db_func)
db_pos = find(p(:,1)==0|p(:,2)==1);
A(db_pos,:) = 0;
A(db_pos,db_pos) = eye(size(db_pos,1));
phi(db_pos) = db_func(p(db_pos,1),p(db_pos,2));
end

