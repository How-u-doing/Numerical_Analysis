% $u = cos(\pi x)\,cos(\pi y)$
% $-\Delta{u} = 2\pi^2 cos(\pi x)\,cos(\pi y)$, in $\Omega$
% $\nabla{u}\cdot\boldsymbol{n}=\pi\,sin(\pi x)\,cos(\pi y)$, on $\Gamma$
% $u=cos(\pi x)\,cos(\pi y)$, on $\partial{\Omega}\ {\backslash}\ {\Gamma}$
% where $\Omega=\{(x,y)\ | \ -\frac{1}{2} < x < 1,\ -1 < y < 1 \}$,
% $\Gamma = \{(x,y)\ | \ x=-\frac{1}{2},\ -1 \leq y \leq 1 \}$

N=6;
L_inf_err=zeros(N,1);
L2_err=zeros(N,1);
H1_semi_err=zeros(N,1);
H1_err=zeros(N,1);
dofs=zeros(N,1);
n=4;
for i=1:N
    [dofs(i),L_inf_err(i),L2_err(i),H1_semi_err(i),H1_err(i)]=run_example1(n);
    n=n*2;
end
L_inf_order=zeros(N-1,1);
L2_order=zeros(N-1,1);
H1_semi_order=zeros(N-1,1);
H1_order=zeros(N-1,1);
for i=1:N-1
    L_inf_order(i)=log2(L_inf_err(i)/L_inf_err(i+1));
    L2_order(i)=log2(L2_err(i)/L2_err(i+1));
    H1_semi_order(i)=log2(H1_semi_err(i)/H1_semi_err(i+1));
    H1_order(i)=log2(H1_err(i)/H1_err(i+1));
end
dofs
disp('L_inf_err   L2_err   H1_semi_err   H1_err')
disp([L_inf_err,L2_err,H1_semi_err,H1_err])
disp('L_inf_order   L2_order   H1_semi_order   H1_order')
disp([L_inf_order,L2_order,H1_semi_order,H1_order])

% plot
figure(3)
% $L^2/L^{\infty}$ norm: $||u-u_h||_*=O(M^{-2})=O(h_M^2)=O(N^{-1})$
% $H^1$ semi-norm: $|u-u_h|_{H^1{(\Omega)}}=O(M^{-1})=O(h_M)=O(N^{-\frac{1}{2}})$
loglog(dofs,L_inf_err,'-^', ...
    dofs,L2_err,'-d', ...
    dofs,exp(-1*log(dofs)+0.8),'--', ...
    dofs,H1_semi_err,'-p', ... dofs,H1_err,'-h', ...
    dofs,exp(-1/2*log(dofs)+1.9),'--');
xlabel('dofs [log]'); xlim([10 3*10^4]);
ylabel('error [log]');
legend('maximum norm','L^2 norm','rate 2','H^1 semi-norm', ... 'H^1 norm',
    'rate 1','Location','SouthWest');
% title('Convergence Rates');
% Since L^2 norm errors are rather small when compared to H^1 semi-norm,
% H^1 semi-norm errors are almost equal to H^1 semi norm, so in the plot
% H^1 semi-norm and H^1 norm almost coincide if we draw them together.

%=========================================================================

function [dofs,L_inf_err,L2_err,H1_semi_err,H1_err] = run_example1(n)
x=[-1/2 1]; y=[-1 1];
% n_x=50; n_y=50;
n_x=n; n_y=n;

[p,t] = generateMesh(x,y,n_x,n_y);
dofs = size(p,1);

% alpha = @(x,y) 1; % alpha = @(x,y) 1*[1 0;0 1];
% gamma = @(x,y) x.*0+y.*0;
A = assembleMatrix(p,t); % A = assembleReactionDiffusionMatrix(p,t,alpha,gamma);
f = @(x,y) 2*pi^2*cos(pi*x).*cos(pi*y);
phi = assembleVector(p,t,f); % phi = assembleVectorByGaussQuad(p,t,f);

% process boundary conditions
nb_func = @(y) -pi*cos(pi*y);
phi = processNeumannBoundary(phi,p,nb_func);
db_func = @(x,y) cos(pi*x).*cos(pi*y);
[A,phi] = processDirichletBoundary(A,phi,p,db_func);

% solve the linear system of equations
u = A\phi;

u_r = @(x, y) cos(pi*x).*cos(pi*y);
grad_u_r = @(x, y) [-pi*sin(pi*x).*cos(pi*y),-pi*cos(pi*x).*sin(pi*y)];
[L_inf_err,L2_err,H1_semi_err,H1_err] = errorEstimates(p,t,u,u_r,grad_u_r);

% t is specified as a 4-by-Nt matrix in pdemesh, although we don't
% use the 4th row in a 2D mesh. To be compatible with 3D meshes?
t=[t';zeros(1,size(t,1))];
p=p'; % required to be 2-by-nP

figure(1)
pdemesh(p,0123,t,u) % 2nd arg `e` is unused, but has to be a matrix
title('LFEM solution')

figure(2)
pdemesh(p,0123,t,u_r(p(1,:),p(2,:)))
title('exact solution')

end

function phi = processNeumannBoundary(phi,p,nb_func)
% $g=\nabla{u}\cdot\boldsymbol{n}=\pi\,sin(\pi x)\,cos(\pi y)$
% $\int_{\Gamma}\,gv\,ds=\sum_i\int_{{\Gamma}_i}\,g\,(v_h^{s(i)}+v_h^{e(i)})\,ds$

pos = find(p(:,1)==-1/2);
edges = [pos(1:end-1),pos(2:end)];
nE = size(edges,1);
for i = 1:nE
    idx_s = edges(i,1); idx_e = edges(i,2);
    y_s = p(idx_s,2); y_e = p(idx_e,2);
    diff_y = y_e-y_s;
    phi(idx_s) = phi(idx_s)+ ...
        diff_y * Gaussquad(@(t)nb_func(y_s+t*diff_y).*t, 0, 1);
    phi(idx_e) = phi(idx_e)+ ...
        diff_y * Gaussquad(@(t)nb_func(y_s+t*diff_y).*(1-t), 0, 1); 
end
end

function [A,phi] = processDirichletBoundary(A,phi,p,db_func)
db_pos = find(p(:,1)==1|p(:,2)==-1|p(:,2)==1);
A(db_pos,:) = 0;
A(db_pos,db_pos) = eye(size(db_pos,1));
phi(db_pos) = db_func(p(db_pos,1),p(db_pos,2));
end

