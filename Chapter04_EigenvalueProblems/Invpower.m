% Determine the eigenvalue of A that is closest to a speci?ed number q. 
% 10170437 Mark Taylor

function [u,x,k]=Invpower(A,x0,tol,N)
% u: (nonzero) dominant eigenval of A,
% x: eigenvector corresponding to u, which satisfies norm(x,inf)=1.

[m,n] = size(A);
if m ~= n 
    error('A must be square!')
end

% Set default function input arguments
if nargin<4
    N=100;
    if nargin<3
        tol=1e-6;
        if nargin<2
            x0=ones(n,1);
            if nargin<1
                error('Too few input argument(s)')
            end            
        end
    end
end

x=x0;
k=1;
% Rayleigh quotient
q=x.'*A*x/(x.'*x);

% Render norm(x,inf)=1.
x=x/x(maxIndex(x,1,n));

u0=0;
u1=0;
while k<=N        
    % Solve the linear system (A?qI)y=x.
    y=GauEli(A-q*eye(n),x);
    
    if ischar(y)==true
        fprintf('q is an eigenvalue');
        u=q;
        return;
    end
    
    % Set u maximum component of y, the dominant eigenvalue converges at u.
    u=y(maxIndex(y,1,n));
    
    % Accelerate convergence via Aitken's delta^2 method 
    u_refine=u0-(u1-u0)^2/(u-2*u1+u0);
    
    err=max(abs(x-y/u));
    x=y/u;
    if abs(err)<tol && k>=4        
        u=u_refine;
        u=1/u+q;
        return;
    end    
    
    k=k+1; 
    u0=u1;
    u1=u;
end

% The number of iterations was exceeded.
k=k-1;% k=N
u=u_refine;
fprintf('\nCannot compute the approximate eigenvalue within %d iterations in the tolerance of %d!\n',N,tol);
fprintf('Calculated eigenvalue & eigenvector in the last iteration are as followed:\n');
    
end
