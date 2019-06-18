% Compute dominant eigenvalue via power method
% 10170437 Mark Taylor

function [u,x,k]=Power(A,x0,tol,N)
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
% Render norm(x,inf)=1.    
x=x/x(maxIndex(x,1,n));    
k=1;
while k<=N
    y=A*x;
    
    % Set u maximum component of y, the dominant eigenvalue converges at u.
    u=y(maxIndex(y,1,n));
    if abs(u)<eps
        fprintf('This matrix has the eigenvalue 0, and its corresponding eigenvector is as followed:\n');
        x        
        prompt='Please select a NEW vector x and restart(by default we provided ones(n,1), please enter a DIFFERENT one)\nx=';
        x=input(prompt);
        while iscolumn(x)==false || length(x)~=n || isequal(x,x0)==true
            fprintf('Invalid x! x must a %d-by-1 vector that differs from x0\n',n);
            x=input(prompt);
        end
        x0=x;        
    end      
    
    err=max(abs(x-y/u));
    x=y/u;
    if abs(err)<tol        
        return;
    end    
    
    k=k+1; 
end

% The number of iterations was exceeded.
k=k-1;% k=N
fprintf('\nCannot compute the approximate dominant eigenvalue within %d iterations in the tolerance of %d!\n',N,tol);
fprintf('Calculated eigenvalue & eigenvector in the last iteration are as followed:\n');
    
end
