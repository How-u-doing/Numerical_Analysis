% Newton Interpolation
% 10170437 Mark Taylor
function s = NI(x,y,xq)
% x,y: vector x & y are sample points and their corresponding function values, respectively
% xq: query points, specified as a vector, to be evaluated
% s: the function values at those interpolating£¨query£© points
n=length(x);
A=NaN(n);   % coefficient matrix of divided differences
A(:,1)=y;
s=y(1);
for j=2:n
    for i=j:n
        A(i,j)=(A(i,j-1)-A(i-1,j-1))/(x(i)-x(i-j+1));
    end
    t=1;
    for k=1:j-1
        t=t.*(xq-x(k));
    end
    s=s+A(j,j).*t;
end
end

