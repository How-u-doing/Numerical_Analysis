% ODE test of several methods
% Example 2.1, P207
f=@(x,y)y-2*x/y;
ya=1;   % Solution: y=sqrt(1+2*t);
a=0;
b=1;
n=10;
[t1, y1]=Forward_Euler(f,a,b,ya,n);
[t2, y2]=Backward_Euler(f,a,b,ya,n);
[t3, y3]=Modified_Euler(f,a,b,ya,n);
[t4, y4]=Midpoint(f,a,b,ya,n);
[t5, y5]=Runge_Kutta(f,a,b,ya,n);

t=0:0.05:1;
y=sqrt(1+2*t);
plot(t,y,'g',t1,y1,'b--',t2,y2,'r--',t3,y3,'c--',t4,y4,'m--',t5,y5,'k--')
legend('Original Function','Forward Euler','Backward Euler','Modified Euler','Midpoint','Runge-Kutta')

figure(2)
[t6, y6]=Trapezoidal(f,a,b,ya,n);
plot(t,y,'g',t3,y3,'b--',t4,y4,'r--',t6,y6,'m--')
legend('Original Function','Modified Euler','Midpoint','Trapezoidal')

