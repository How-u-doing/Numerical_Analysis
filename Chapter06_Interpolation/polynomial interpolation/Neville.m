% Neville's Method
% 10170437 Mark Taylor
function [s, Q] = Neville(x,y,xq)
% x,y: vector x & y are sample points and their corresponding function values, respectively
% xq: query points, specified as a vector, to be evaluated
% s: the function values at those interpolating£¨query£© points
n=length(x);
m=length(xq);
s=zeros(1,m);
Q=NaN(n);
Q(:,1)=y;
for k=1:m
    for j=2:n
        for i=j:n
            Q(i,j)=((xq(k)-x(i-j+1))*Q(i,j-1)-(xq(k)-x(i))*Q(i-1,j-1))/(x(i)-x(i-j+1));
        end    
    end
    s(k)=Q(n,n);
end

end