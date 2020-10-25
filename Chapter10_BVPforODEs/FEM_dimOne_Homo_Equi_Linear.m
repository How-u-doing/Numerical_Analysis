function [A,F,u] = FEM_dimOne_Homo_Equi_Linear(p,q,f,a,b,n)
%=================================================================
% One dimensional two - point Homogeneous boundary value problem
% Equidistant subdivision
% FEM_Piecewise linear interpolation
% -d(p * du/dx)/dx + q*u = f
%  u(a) = 0;  u'(b) = 0
% ***********************************
%   [a,b]上n等分,分段线性插值
% h: Equidistant step
% x: nodes
% ***********************************
%   A*(u0;u1;...;un) = (f0;f1;...;fn)
%   u0=0
%  故求解(u1;...;un) 按A(2:n,2:n)\（f1;...;fn)
%======================================================

h = (b-a)/n;
x = a:h:b;
A = zeros(n+1); % 计算时第一行第一列将删掉
F = zeros(n+1,1);

for i = 1:n
    A1 = @(k) p(x(i)+h*k)/h + h*q(x(i)+h*k).*(1-k).^2;
    A2 = @(k) p(x(i)+h*k)/h + h*q(x(i)+h*k).*(k).^2;
    A3 = @(k) -p(x(i)+h*k)/h + h*q(x(i)+h*k).*k.*(1-k);
    F1 = @(k) f(x(i)+h*k).*(1-k);
    F2 = @(k) f(x(i)+h*k).*k;
    IA1 = integral( A1 , 0 , 1); % A(i,i)_i
    IA2 = integral( A2 , 0 , 1); % A(i+1,i+1)_i
    IA3 = integral( A3 , 0 , 1); % A(i+1,i)_i = a(i,i+1)_i
    IF1 = h * integral( F1 , 0 , 1);
    IF2 = h * integral( F2 ,0 , 1);

    A(i,i) = A(i,i)+IA1;
    A(i+1,i+1) = A(i+1,i+1)+IA2;
    A(i,i+1) = A(i,i+1)+IA3;
    A(i+1,i) = A(i+1,i)+IA3;
    F(i) = F(i)+IF1;
    F(i+1) = F(i+1)+IF2;
    
end
A = A(2:n+1,2:n+1);
F = F(2:n+1);
u = A\F;

end