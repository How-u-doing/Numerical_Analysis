% Homework( Exercise 1 ), P67
% 10170437 Mark Taylor
x=-1:2;
y=[4 1 6 26];
dya=-6;dyb=15;
s=CShw(x,y,dya,dyb);
syms t;
X=zeros(3,11);Y=X;
for i=1:3
    fprintf('%d<=t<=%d:\n',i-2,i-1);    
    S=s(i)
    X(i,:)=i-2:.1:i-1;
    Y(i,:)=subs(S,t,X(i,:));
end
X=[X(1,:),X(2,:),X(3,:)];
Y=[Y(1,:),Y(2,:),Y(3,:)];
plot(x,y,'o',X,Y,'r--')
legend('Sample points','Cubic Spline')


