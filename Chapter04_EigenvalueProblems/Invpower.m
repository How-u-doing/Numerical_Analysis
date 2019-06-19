% Determine the eigenvalue of A that is closest to a specified number q,
% where q is an approximate eigenvalue of A, and we strongly suggest you 
% set it 'default' or 'max' if you wanna get the spectral radius of A.
% This method is often used to refine an approximate eigenvalue and to compute
% an eigenvector once an eigenvalue has been found by some other technique. 
% 10170437 Mark Taylor

function [u,x,k]=Invpower(A,q,x,tol,N)
% u: eigenval of A that is "closest" to the specified number q,
% x: eigenvector corresponding to u, which satisfies norm(x,inf)=1.


% /*Something worthing noting: 
%
% u is definitely one of eigenvals of A, but is not absolutely the closest 
% to q, it depends on tol. The smaller tol is, the closer u is approaching
% q. And if q is about in the meddle of two eigenvalues of A, the iterative 
% times is very large! See "test1.m" to verify the fact.
%
% **/


[m,n] = size(A);
if m ~= n 
    error('A must be square!')
end

% Use flag that tells if to solve the maximum eigenvalue in absolute value
% to quickly solve the objective eigenvalue & eigenvector. 
isSolvingMax=false;

% Set default function input arguments
if nargin<5
    N=100;
    if nargin<4
        tol=1e-6;
        if nargin<3
            x=ones(n,1);
            if nargin<2
                % We set q Rayleigh quotient by default.
                % q->u_max(in absolute value)
                q=x.'*A*x/(x.'*x);
                isSolvingMax=true;
                if nargin<1
                    error('Too few input argument(s)!')
                end
            end            
        end
    end
end


if ischar(q)
    if strcmp(q,'default')==true || strcmp(q,'max')==true
        q=x.'*A*x/(x.'*x);
        isSolvingMax=true;
    else
        error('Invalid q! q must be numeric or ''default'' or ''max''.')
    end
end

% Render norm(x,inf)=1.
x=x/x(maxIndex(x,1,n));

k=1;
while k<=N        
    % Solve the linear system (A-qI)y=x.
    y=GauEli(A-q*eye(n),x);
    
    if ischar(y)
        fprintf('\nq is an eigenvalue.\n');
        u=q;
        return;
    end
    
    % Set u maximum component of y, the dominant eigenvalue converges at u.
    u=y(maxIndex(y,1,n));
    
    
    err=max(abs(x-y/u));
    x=y/u;
    
    if abs(err)<tol        
        u=1/u+q;
        return;
    end 
    
    % Accelerate convergence via Rayleigh quotient when we 
    % wanna solve the maximum eigenvalue in absolute value.
    if isSolvingMax==true
        q=x.'*A*x/(x.'*x);
    end
    
    k=k+1;
end

% The number of iterations was exceeded.
k=k-1;% k=N
u=1/u+q;
fprintf('\nCannot compute the approximate eigenvalue within %d iterations in the tolerance of %d!\n',N,tol);
fprintf('Calculated eigenvalue & eigenvector in the last iteration are as followed:\n');
    
end
