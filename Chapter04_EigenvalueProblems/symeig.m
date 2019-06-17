% Compute m-th largest eigenvalue of symmetric matrices via bisection method.
% In our text, we call it Givens-Householder method.  
% 10170437 Mark Taylor

function [eigval,k]=symeig(A,m,a,b,tol,N)

if isequal(A,A.')==false 
    error('A must be symmetric!')
end

% Upper-Hessenberg reduction
C=Hessenberg(A); % C is a tridiagonal matrix who has identical eigenvalues of A

n = size(C,2);
% Set default function input arguments
if nargin<6
    N=100;
    if nargin<5
        tol=1e-6;
        if nargin<4
            b=norm(C,inf);
            if nargin<3
                a=-b;
                if nargin<2
                    m=1; % Solve the maximum eigenvalue by default
                    if nargin<1
                        error('Too few input argument(s)')
                    end
                end
            end
        end
    end
end


if strcmp(b,'default')==true
    b=norm(C,inf);
end
if strcmp(a,'default')==true
    a=-b;
end
if strcmp(m,'max')==true
    m=1;
end
if strcmp(m,'min')==true
    m=n;
end

% Check the validity of input arguments
if a>=b
    error('a must be less than b, i.e. a<b')    
end
if m<1 || m>n
    error('m must satisfy 1<=m<=%d',n)    
end



[eigval,k]=bisection(C,m,a,b,tol,N);

end