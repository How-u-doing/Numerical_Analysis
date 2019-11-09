% Exer 11, page 138 of our text
x=[1 3 4 6 7];
y=[-2.1 -0.9 -0.6 0.6 0.9];
rho=[0.1 0.2 0.4 0.1 0.2];
N=length(x);
n=2;              % approximation degree
a=lsq(x,y,n,rho);
% a=lsq0(x,y,n,rho);
p=flip(a);
xx=1:.01:7;
yls=polyval(p,xx);
plot(x,y,'o',xx,yls,'r--');
xlabel('x');ylabel('y');
legend('Captured Data','Least Squares')


