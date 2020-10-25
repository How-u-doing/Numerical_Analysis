% Piecewise Quadratic element Rayleigh-Ritz method
% element stiffness matrix version
% Inspired by Miss White ;)
function u = PQRR_es(f,p,q,a,b,n)

h = (b-a)/n; % step length. no check for N, hhh
A = zeros(2*n+1);   % coef matrix, symmetric
F = zeros(2*n+1,1); % right-hand side (RHS) vector
x = linspace(a,b,n+1);

% k = ( x- x_{j-1} ) / h_j, 0 <= k <=1
% u_h = u_{j-1}*(2k-1)(k-1) + u_{j-1/2}*4k(1-k) + u_j*k(2k-1)
% u_h' = u_{j-1}*(4k-3)/h_j + u_{j-1/2}*(-8k+4)/h_j + u_j*(4k-1)/h_j

% hopefully inlined
phi1 = @(k)(2*k-1).*(k-1);
phi2 = @(k)4*k.*(1-k);
phi3 = @(k)k.*(2*k-1);
k1 = @(k)(4*k-3)/h;
k2 = @(k)(4-8*k)/h;
k3 = @(k)(4*k-1)/h;

K = zeros(3); % element stiffness matrix
for j = 1 : n
    K(1,1) = integral(@(k)(p(k*h+x(j)).*k1(k).^2+...
        q(k*h+x(j)).*phi1(k).^2)*h, 0, 1); % omit *h in each item for optimization
    
    K(1,2) = integral(@(k)(p(k*h+x(j)).*k1(k).*k2(k)+...
        q(k*h+x(j)).*phi1(k).*phi2(k))*h, 0, 1);
    
    K(1,3) = integral(@(k)(p(k*h+x(j)).*k1(k).*k3(k)+...
        q(k*h+x(j)).*phi1(k).*phi3(k))*h, 0, 1);
    
    K(2,2) = integral(@(k)(p(k*h+x(j)).*k2(k).^2+...
        q(k*h+x(j)).*phi2(k).^2)*h, 0, 1);
    
    K(2,3) = integral(@(k)(p(k*h+x(j)).*k2(k).*k3(k)+...
        q(k*h+x(j)).*phi2(k).*phi3(k))*h, 0, 1);
    
    K(3,3) = integral(@(k)(p(k*h+x(j)).*k3(k).^2+...
        q(k*h+x(j)).*phi3(k).^2)*h, 0, 1);
    
    % using symmetricity
    K(2,1) = K(1,2);
    K(3,1) = K(1,3);
    K(3,2) = K(2,3);
    
    % assemble the element stiffness matrix
    A(2*j-1,2*j-1) = A(2*j-1,2*j-1) + K(1,1);
    A(2*j-1, 2*j ) = A(2*j-1, 2*j ) + K(1,2);
    A(2*j-1,2*j+1) = A(2*j-1,2*j+1) + K(1,3);
    A( 2*j ,2*j-1) = A( 2*j ,2*j-1) + K(2,1);
    A( 2*j , 2*j ) = A( 2*j , 2*j ) + K(2,2);
    A( 2*j ,2*j+1) = A( 2*j ,2*j+1) + K(2,3);
    A(2*j+1,2*j-1) = A(2*j+1,2*j-1) + K(3,1);
    A(2*j+1, 2*j ) = A(2*j+1, 2*j ) + K(3,2);
    A(2*j+1,2*j+1) = A(2*j+1,2*j+1) + K(3,3);
    
    % calculate & (implicitly) assemble right-hand side vector
    F(2*j-1) = F(2*j-1) + integral(@(k)(f(k*h+x(j)).*phi1(k))*h, 0, 1);
    F( 2*j ) = F( 2*j ) + integral(@(k)(f(k*h+x(j)).*phi2(k))*h, 0, 1);
    F(2*j+1) = F(2*j+1) + integral(@(k)(f(k*h+x(j)).*phi3(k))*h, 0, 1);
end

A = A(2:2*n+1,2:2*n+1);
F = F(2:2*n+1);

u = A\F;
end


