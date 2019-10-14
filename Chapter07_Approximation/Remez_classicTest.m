% Remez's 1st algorithm test on 1/(1+x^2)
% 10170437 Mark Taylor
f=@(x) 1./(1+x.^2);
a=-2;b=2;
syms t;
% click <Run> or umcomment others to focus on a particular test
% you can also just choose one test and change f,a,b,x to play around

% % linear approximation test
% x1=[-1,0,1];
% y1=f(x1);
% [pn1,Enf1,x1,k1]=Remez(f,x1,y1,a,b)
% % you can verify the validity like this
% g1=pn1-f(t);
% vpa(subs(g1,t,x1),5)
% xx=a:.02:b;
% ypn1=subs(pn1,t,xx);
% figure
% plot(xx,f(xx),'g',xx,ypn1,'--r')
% legend('Original Function','Optimal Polynomial')
% title('Linear Optimal Polynomial')

% % quadratic approximation test
% x2=[0.1,0.3,0.7,0.9];
% y2=f(x2);
% [pn2,Enf2,x2,k2]=Remez(f,x2,y2,a,b) 
% g2=pn2-f(t);
% vpa(subs(g2,t,x2),5)
% xx=a:.02:b;
% ypn2=subs(pn2,t,xx);
% figure
% plot(xx,f(xx),'g',xx,ypn2,'--r')
% legend('Original Function','Optimal Polynomial')
% title('Quadratic Optimal Polynomial')

% higher degree (7th degree) optimal polynomial test. It needs some time!
x7=[-2,-1,-0.8,-0.5,0,0.3,0.8,1.3,1.7];
y7=f(x7);
xx=a:.02:b;
tt=linspace(-2,2,7);
ynew=NI(tt,f(tt),xx); % for Newton's 
[pn7,Enf5,x7,k7]=Remez(f,x7,y7,a,b)
g7=pn7-f(t);
vpa(subs(g7,t,x7),5)
ypn7=subs(pn7,t,xx);
figure
plot(xx,f(xx),'g',xx,ypn7,'--r',xx,ynew,'--b')
legend('Original Function','Optimal Polynomial','Newton''s P7')
title('7th-degree Optimal Polynomial')


