% Explicit Midpoint method (Runge-Kutta method of order 2)
function [t, y]=Midpoint(f,a,b,ya,n)
h=(b-a)/n;
t=a:h:b;
y=zeros(1,n+1);
y(1)=ya;
for i=2:n+1
    y(i)=y(i-1)+h*f(t(i-1)+h/2,y(i-1)+h/2*f(t(i-1),y(i-1)));    
end
end

