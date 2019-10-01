% B-Spline base function vector form test
x=-5:5;
k=3;
t=-5:0.01:5;
n=length(x);
y=zeros(1,n);
for i=1:n-k-1
    y = BSbaseVec(i,k,x,t);
    plot(t,y);
    pause(0.5);
end