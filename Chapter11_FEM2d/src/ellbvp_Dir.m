% FEM for Elliptic BVP
% $-\nabla\cdot(\bm{\alpha}\nabla u) + \gamma u = f$, in $\Omega$,
% $u = g_d$, on $\partial\Omega$,
% where $\Omega=(0,1)\times(0,1)$.

u_case = 1; % <--- choose which one to run
N=6;
n=4; % $h(\mathcal{M}):=n^{-1}$
n_arr = 2.^(2:N+1).';
L_inf_err=zeros(N,1);
L2_err=zeros(N,1);
H1_semi_err=zeros(N,1);
H1_err=zeros(N,1);
for i=1:N
    [~,L_inf_err(i),L2_err(i),H1_semi_err(i),H1_err(i)]=run_ellbvp(n,u_case);
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
p_L_inf = polyfit(log2(n_arr),log2(L_inf_err),1);
p_L2 = polyfit(log2(n_arr),log2(L2_err),1);
p_H1_semi = polyfit(log2(n_arr),log2(H1_semi_err),1);
p_H1 = polyfit(log2(n_arr),log2(H1_err),1);
n_arr
disp('L_inf_err   L2_err   H1_semi_err   H1_err')
disp([L_inf_err,L2_err,H1_semi_err,H1_err])
disp('L_inf_order   L2_order   H1_semi_order   H1_order')
disp([L_inf_order,L2_order,H1_semi_order,H1_order])
disp('Orders by least square approximation')
disp(-[p_L_inf(1),p_L2(1),p_H1_semi(1),p_H1(1)])

% plot
figure(3)
b = [1.9,0.8; -0.7,-1.7; 5.5,4.3];
% $L^2/L^{\infty}$ norm: $||u-u_h||_*=O(h_\mathcal{M}^2)=O(n^{-2})$
% $H^1$ semi-norm: $|u-u_h|_{H^1{(\Omega)}}=O(h_\mathcal{M})=O(n^{-1})$
loglog(n_arr,L_inf_err,'-^', ...
    n_arr,L2_err,'-d', ...
    n_arr,exp(-2*log(n_arr)+b(u_case,2)),'--', ... % adjust the coef b in exp to
    n_arr,H1_err,'-p', ...                         % make the ref line work better
    n_arr,exp(-1*log(n_arr)+b(u_case,1)),'--');
xlabel('n [log]');
ylabel('error [log]');
legend('maximum norm','L^2 norm','rate 2','H^1 norm', ...
    'rate 1','Location','SouthWest');
% title('Convergence Rates');

%=========================================================================

function [dofs,L_inf_err,L2_err,H1_semi_err,H1_err] = run_ellbvp(n, u_case)
x=[0 1]; y=[0 1];
% n_x=50; n_y=50;
n_x=n; n_y=n;

[p,t] = generateMesh(x,y,n_x,n_y);
dofs = size(p,1);

if u_case == 1
    u_r = @(x, y) sin(pi*x).*cos(pi*y);
    grad_u_r = @(x, y) [pi*cos(pi*x).*cos(pi*y), -pi*sin(pi*x).*sin(pi*y)];    
    alpha = @(x,y)[x.^2+y.^2+1,x.*y; x.*y,x.^2+y.^2+1];
    gamma = @(x,y) 0.*x.*y;
    f = @(x,y) -pi*( 3*x.*cos(pi*x).*cos(pi*y)-2*pi*x.*y.*cos(pi*x).*sin(pi*y)- ...    
               2*pi*(x.^2+y.^2+1).*sin(pi*x).*cos(pi*y)-3*y.*sin(pi*x).*sin(pi*y) );
%     alpha = @(x,y) 1; % alpha = @(x,y) 1*[1 0;0 1];
%     gamma = @(x,y) x.*0+y.*0;
%     f = @(x,y) 2*pi^2*sin(pi*x).*cos(pi*y);
elseif u_case == 2
    u_r = @(x, y) x.*(1-x).*y.*(1-y);
    grad_u_r = @(x, y) [(1-2*x).*y.*(1-y), x.*(1-x).*(1-2*y)];
    % Here the diffusion coefficient $\bm{\alpha}$ is not uniformly positive
    % definite, so the rates of convergence (in $L^\infty, L^2, H^1$ norms)
    % are unspecified though the exact solution is smooth.
    alpha = @(x,y) x.*y; % alpha = @(x,y) x.*y.*[1 0;0 1];
    gamma = @(x,y) 0.*x.*y;
    f = @(x,y) -( (1-4*x).*y.^2.*(1-y)+x.^2.*(1-x).*(1-4*y) );
elseif u_case == 3
    k = 3; % k = 9;
    u_r = @(x, y) sin(2*pi*x).*sin(2*k*pi*y);
    grad_u_r = @(x, y) [2*pi*cos(2*pi*x).*sin(2*k*pi*y), 2*k*pi*sin(2*pi*x).*cos(2*k*pi*y)];
    alpha = @(x,y) [k^2, 0; 0, 1];
    gamma = @(x,y) 0.*x.*y;
    f = @(x,y) 8*k^2*pi^2*sin(2*pi*x).*sin(2*k*pi*y);
end

A = assembleReactionDiffusionMatrix(p,t,alpha,gamma);
% phi = assembleVector(p,t,f);
phi = assembleVectorByGaussQuad(p,t,f);

% process boundary conditions
% see also "example1.m" for how to process Neumann BCs
[A,phi] = processDirichletBoundary(A,phi,p,u_case);

% solve the linear system of equations
u = A\phi;

% get discretization errors in various norms
[L_inf_err,L2_err,H1_semi_err,H1_err] = errorEstimates(p,t,u,u_r,grad_u_r);

% t is specified as a 4-by-Nt matrix in pdemesh, although we don't
% use the 4th row in a 2D mesh. To be compatible with 3D meshes?
t=[t';zeros(1,size(t,1))];
p=p'; % required to be 2-by-nP

figure(1)
pdemesh(p,0123,t,u) % 2nd arg `e` is unused, but has to be a matrix
xlabel('x'); ylabel('y');
title('LFEM solution')

figure(2)
pdemesh(p,0123,t,u_r(p(1,:),p(2,:)))
xlabel('x'); ylabel('y');
title('exact solution')
end

function [A,phi] = processDirichletBoundary(A,phi,p,u_case)
if u_case == 1
    db_pos1 = find(p(:,1)==0|p(:,1)==1);
    db_pos2 = find(p(:,2)==0);
    db_pos3 = find(p(:,2)==1);
    db_pos = [db_pos1; db_pos2; db_pos3];
    db_pos = unique(db_pos);
    phi(db_pos1) = zeros(size(db_pos1,1),1);
    phi(db_pos2) = sin(pi*p(db_pos2,1));
    phi(db_pos3) = -sin(pi*p(db_pos3,1));
else%if u_case == 2 || u_case ==3
    db_pos = find(p(:,1)==0|p(:,1)==1|p(:,2)==0|p(:,2)==1);    
    phi(db_pos) = 0;
end
A(db_pos,:) = 0;
A(db_pos,db_pos) = eye(size(db_pos,1));
end


