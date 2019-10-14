% Remez's 1st algorithm test on polynomials
% Of course, when no. of test points is no less than the degree of f(x) plus
% 2, we can get goddamn perfect results, but we can use other more effective
% algorithms to determine such pricise polynomial. Our interest is, however, 
% try to find a lower degree optimal polyniomial to approximate a higher one.
% 10170437 Mark Taylor
f=@(x)x.^5-5*x.^2+9;% Try: change to 2*x.^3+x.^5-x.^6-x+6, 3*x.^2./(5-x)
a=-2;b=2;
syms t;
% click <Run> or umcomment others to focus on a particular test
% you can also just choose one test and change f,a,b,x to play around

% linear approximation test
x1=[-1,0,1];
y1=f(x1);
[pn1,Enf1,x1,k1]=Remez(f,x1,y1,a,b)
% you can verify the validity like this
g1=pn1-f(t);
vpa(subs(g1,t,x1),5)
xx=a:.02:b;
ypn1=subs(pn1,t,xx);
figure
% When f''(x)>0 (e.g. x.^4+x.^2-7), -1->a, 1->b, and '--r' is parallel to line 'b'
plot(xx,f(xx),'g',xx,ypn1,'--r',[a b],[f(a),f(b)],'b')
legend('Original Function','Optimal Polynomial')
title('Linear Optimal Polynomial')

% quadratic approximation test
x2=[0.1,0.3,0.7,0.9];
y2=f(x2);
[pn2,Enf2,x2,k2]=Remez(f,x2,y2,a,b) 
g2=pn2-f(t);
vpa(subs(g2,t,x2),5)
xx=a:.02:b;
ypn2=subs(pn2,t,xx);
figure
plot(xx,f(xx),'g',xx,ypn2,'--r')
legend('Original Function','Optimal Polynomial')
title('Quadratic Optimal Polynomial')

% higher degree (5th degree) optimal polynomial test
x5=[-1,-0.5,0,0.3,0.4,0.6,0.7];
y5=f(x5);
[pn5,Enf5,x5,k5]=Remez(f,x5,y5,a,b)
g5=pn5-f(t);
vpa(subs(g5,t,x5),5)
xx=a:.02:b;
ypn5=subs(pn5,t,xx);
figure
plot(xx,f(xx),'g',xx,ypn5,'--r')
legend('Original Function','Optimal Polynomial')
title('5th-degree Optimal Polynomial')


