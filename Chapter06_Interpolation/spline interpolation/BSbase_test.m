% B-Spline base function(component-wise) test
x=-5:5;
k=3;
t=-5:0.01:5;
n=length(x);
N=length(t);
y=zeros(1,N);
for i=1:n-k-1
    for j=1:N
        y(j) = BSbase(i,k,x,t(j));
    end
    hold on
    plot(t,y);
    pause(0.5);
end
