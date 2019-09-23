% Homework( Exercise 1 ), P67
% 10170437 Mark Taylor
x=-1:2;
y=[4 1 6 26];
dya=-6;dyb=15;
s=CShw(x,y,dya,dyb);
for i=1:3
    fprintf('%d<=t<=%d:\n',i-2,i-1);
    S=s(i)
end

