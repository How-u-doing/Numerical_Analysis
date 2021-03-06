% Romberg Integration (iterative method)
function [T,k]=Romberg(f,a,b,maxIt,tol)
% set default args
if nargin<5
    tol=1e-6;
    if nargin<4
        maxIt=30;   % 2^30 intervals at most
        if nargin<3
            error('Too few arguments')
        end
    end
end

h=(b-a);
% R1, R2 are used to store the most recent rows
R1=zeros(1,maxIt+1);
R1(1)=h/2*(f(a)+f(b));
for k=1:maxIt
    n=2^(k-1);
    h=(b-a)/n; % h is the even interval of R1
    H=0;
    for i=0:n-1
        % when n is large, say 2^25, it's computation expensive 
        H=H+f(a+(i+1/2)*h);
    end
    H=h*H;
    
    R2=zeros(1,maxIt+1);
    R2(1)=1/2*(R1(1)+H);
    for j=2:k+1
        R2(j)=1/(4^(j-1)-1)*(4^(j-1)*R2(j-1)-R1(j-1));
    end
    
    % tell when to stop
    if abs(R2(k+1)-R1(k))<tol
        T=R2(k+1);
        return
    end
    
    % update R1
    R1=R2;
end
% Maximum number of iterations exceeded (the procedure was successful).
fprintf('Maximum no. of iterations exceeded, the results of last iteration are as follows:\n')
T=R2(maxIt+1);
end

