% Remez's 1st algorithm test on sin(x) over [0,pi] 
% 10170437 Mark Taylor
f=@(x)sin(x);
a=0;b=pi;
syms t;
% click <Run> or umcomment others to focus on a particular test
% you can also just choose one test and change f,a,b,x to play around

% linear approximation test
% change x1 to [1,2,3] to see what's going on (1->a, 3->b, since f''(x)<=0) 
x1=[0,2,pi]; % you can even change 2 to 0.01 !
y1=f(x1);
[pn1,Enf1,x1,k1]=Remez(f,x1,y1,a,b)
% you can verify the validity like this
g1=pn1-f(t);
vpa(subs(g1,t,x1),5)
xx=a:.1:b;
ypn1=subs(pn1,t,xx);figure
plot(xx,f(xx),'g',xx,ypn1,'--r')
legend('Original Function','Optimal Polynomial')
title('Linear Optimal Polynomial')

% quadratic approximation test
x2=[0,1,3,pi];
y2=f(x2);
% There would have some warnings when we attempt to solve its sup point by 
% finding all extreme points and comparing those extremum values. Ignore them!
[pn2,Enf2,x2,k2]=Remez(f,x2,y2,a,b) 
g2=pn2-f(t);
vpa(subs(g2,t,x2),5)
xx=a:.1:b;
ypn2=subs(pn2,t,xx);figure
plot(xx,f(xx),'g',xx,ypn2,'--r')
legend('Original Function','Optimal Polynomial')
title('Quadratic Optimal Polynomial')

% higher degree (5th degree) optimal polynomial  test
x5=[0,1,1.5,2,2.5,3,pi];% change 0 to .5, pi to 3.1, see what's going on
y5=f(x5);
[pn5,Enf5,x5,k5]=Remez(f,x5,y5,a,b)
g5=pn5-f(t);
vpa(subs(g5,t,x5),5)
xx=a:.1:b;
ypn5=subs(pn5,t,xx);figure
plot(xx,f(xx),'g',xx,ypn5,'--r')
legend('Original Function','Optimal Polynomial')
title('5th-degree Optimal Polynomial')


