% Bisection Mothod for the eigenvalues of tridiagonal matrix 
% 10170437 Mark Taylor

function [eigval,k]=bisection(C,m,a,b,tol,N)
% Invoke this function only when C is confirmed tridiagonal.
% This funtion is used to solve eigenvalue(s) range from a to b, i.e. [a,b)
% Let u_i be one of C's eigenvalues, i.e. C*x=u_i*x, i=1,2,...n.
% We can use specific range [a,b) to accelerate convergence, but be careful 
% that u_m must be in the range of [a,b), otherwise we only get a or b that 
% is closer to genuine eigenvalue as the (most likely incorrect)answer. If   
% you are not sure about it, we strongly suggest you set a and b 'default'.

% u_1>=...>=u_(m-1)>=u_m>=u_(m+1)>=...>=u_n


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


r=(a+b)/2;
k=1;
while k<=N    
    nNonnegative=0;   
    for i=1:n
        if q(C,r,i)>=0
            nNonnegative=nNonnegative+1;
        end
    end
    
    if nNonnegative>=m
        a=r;
    else
        b=r;
    end
    
    r=(a+b)/2;
    if abs(b-a)<tol
        eigval=r;
        return;
    end
    
    k=k+1; 
end

% The number of iterations was exceeded.
k=k-1;% k=N
eigval=r;
fprintf('\nCannot compute the approximate eigenvalue within %d iterations in the tolerance of %d!\n',N,tol);
fprintf('Calculated eigenvalue in the last iteration is as followed:\n');
end



function y=q(C,u,k)

if k==0
    y=1;
    return;
end

if k==1
    y=C(1,1)-u;
    return;
end

% 2<=k<=n
if q(C,u,k-1)*q(C,u,k-2)~=0
    y=C(k,k)-u-C(k,k-1)^2/q(C,u,k-1);
    return;
end

if q(C,u,k-2)==0
    y=C(k,k)-u;
    return;
end

if q(C,u,k-1)==0
    y=-1;
    return;
end

end

