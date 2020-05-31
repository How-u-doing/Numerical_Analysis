% Runge-Kutta method of order 4
function [t, y]=Runge_Kutta(f,a,b,ya,n)
h=(b-a)/n;
t=a:h:b;
y=zeros(1,n+1);
y(1)=ya;
for i=2:n+1
    k1=f(t(i-1),y(i-1));
    k2=f(t(i-1)+h/2,y(i-1)+h/2*k1);
    k3=f(t(i-1)+h/2,y(i-1)+h/2*k2);
    k4=f(t(i-1)+h,y(i-1)+h*k3);
    y(i)=y(i-1)+h/6*(k1+2*k2+2*k3+k4);
end
end
