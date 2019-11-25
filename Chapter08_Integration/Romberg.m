% Romberg Integration £¨iterative method)
function [T,k]=Romberg(f,a,b,maxIt,tol)
% set default args
if nargin<5
    tol=1e-6;
    if nargin<4
        maxIt=100;
        if nargin<3
            error('Too few arguments')
        end
    end
end

h=(b-a);
% R1, R2 are used to store the most recent rows
R1=zeros(1,maxIt);
R1(1)=h/2*(f(a)+f(b));
R2=R1;
for k=1:maxIt
    n=2^k;
    h=(b-a)/n;
    H=0;
    for i=0:n-1
        H=H+f(a+(i+1/2)*h);
    end
    H=h*H;
    
    R2(1)=1/2*(R1(1)+H);
    for j=2:k+1
        R2(j)=1/(4^k-1)*(4^k*R2(j-1)-R1(j-1));
    end
    
    % tell when to stop
    if abs(R2(1)-R1(1))<tol
        T=R2(k+1);
        return
    end
    
    % update R1
    R1=R2;
end

end

