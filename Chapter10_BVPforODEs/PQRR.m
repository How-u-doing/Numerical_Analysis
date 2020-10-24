% Piecewise Quadratic element Rayleigh-Ritz method
% ===============================================================
% Solve a linear two-point boundary-value problem
%                   -(p(x)u'(x))'+q(x)u(x)=f(x),    a<=x<=b,
% with initial conditions
%                           u(a)=0, u'(b)=0.
% ===============================================================
% See the theory part under folder <quadratic_element_img>
function uc=PQRR(f,p,q,a,b,N)
% INPUT:
%   f: f(x)
%   p: p(x); optional, by default p(x)=1
%   q: q(x); optional, by default q(x)=0
%   a, b: interval [a,b]; optional, by default [a,b]=[0,1]
%   N: # of evenly spaced intervals; optional, by default N=10
% OUTPUT:
%   uc: coef vector of approx solution

% ===============================================================
% *****  NOTE that f, p, q must be element-wise functions  *****
% e.g. f(x)=@(x)x.^2+1 (NOT x^2+1)
% ===============================================================

% set default args
if nargin<6
    N=10;
    if nargin<5
        b=1;
        if nargin<4
            a=0;
            if nargin<3
                q=@(x)x.*0;
                if nargin<2
                    p=@(x)x.^0;
                    if nargin<1                 
                        error('Error! Function f is not set.')
                    end
                end
            end
        end
    end
end
if b<a
    tmp=a; a=b; b=tmp;
end

h=(b-a)/N; % step length. no check for N, hhh
A=zeros(2*N);   % coef matrix, symmetric
c=zeros(2*N,1); % right-hand side (RHS) vector
xx=linspace(a,b,N+1);

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

% determine tridiagonal elements of A, and RHS vector c
for j = 1 : N
    % Or use MATLAB func <integral> instead of <Gaussquad>.
    % --> <Ctrl + F> --> Replace All. Then go to PQRR_test
    % script file & run it again, it shall get same graph.
    A(2*j-1,2*j-1) = Gaussquad(@(k)(p(k*h+xx(j)).*k2(k).^2+...
        q(k*h+xx(j)).*phi2(k).^2)*h, 0, 1); % omit *h in each item for optimization
    A(2*j-1,2*j) = Gaussquad(@(k)(p(k*h+xx(j)).*k2(k).*k3(k)+...
        q(k*h+xx(j)).*phi2(k).*phi3(k))*h, 0, 1);
    A(2*j,2*j-1) =  A(2*j-1, 2*j ); % symmetric
    
    A(2*j,2*j) = Gaussquad(@(k)(p(k*h+xx(j)).*k3(k).^2+...
        q(k*h+xx(j)).*phi3(k).^2)*h, 0, 1);
    if(j<N)
        A(2*j,2*j)=A(2*j,2*j)+Gaussquad(@(k)(p(k*h+xx(j+1)).*k1(k).^2+...
            q(k*h+xx(j+1)).*phi1(k).^2)*h, 0, 1);
        
        A(2*j,2*j+1) = Gaussquad(@(k)(p(k*h+xx(j+1)).*k1(k).*k2(k)+...
            q(k*h+xx(j+1)).*phi1(k).*phi2(k))*h, 0, 1);
        
        A(2*j,2*j+2) = Gaussquad(@(k)(p(k*h+xx(j+1)).*k1(k).*k3(k)+...
            q(k*h+xx(j+1)).*phi1(k).*phi3(k))*h, 0, 1);
        
        A(2*j+1,2*j) = A( 2*j ,2*j+1);
        A(2*j+2,2*j) = A( 2*j ,2*j+2);
    end
    
    c(2*j-1) = Gaussquad(@(k)(f(k*h+xx(j)).*phi2(k))*h, 0, 1);
    c( 2*j ) = Gaussquad(@(k)(f(k*h+xx(j)).*phi3(k))*h, 0, 1);
    if(j<N)
        c(2*j)=c(2*j)+Gaussquad(@(k)(f(k*h+xx(j+1)).*phi1(k))*h, 0, 1);
    end
end

% uc=A\c; % coef vector of approx solution
% or use user-defined function <GauEli>
uc=GauEli(A,c);

end



