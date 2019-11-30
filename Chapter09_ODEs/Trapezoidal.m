% Implicit trapezoidal rule 
function [t, y]=Trapezoidal(f,a,b,ya,n,maxIt,tol)
if nargin<7
    maxIt=1000;
    if nargin<6
        tol=1e-6;
        if nargin<5
            error('Too few arguments')
        end
    end 
end
h=(b-a)/n;
t=a:h:b;
y=zeros(1,n+1);
y(1)=ya;
for i=2:n+1
    y1=y(i-1)+h*f(t(i-1),y(i-1));
    % iterate to solve y(i)
    for k=1:maxIt
        y2=y(i-1)+h/2*(f(t(i-1),y(i-1))+f(t(i),y1));
        if abs(y1-y2)<tol
            y(i)=y2;
            break;
        end
        % update y1
        y1=y2;
    end
end
end

