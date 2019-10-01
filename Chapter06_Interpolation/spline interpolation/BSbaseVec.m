% Attempt to construct a B-spline base function via vector recursion
% 10170437 Mark Taylor
function y = BSbaseVec(i,k,x,t)
% i,k: the B-spline base of kth order(degree of k-1) whose support lies in interval (x_i,x_(i+k))
%  x : sample points, [x_1,x_2,...,x_n]
%  t : points to be evaluated
m=length(t);
y=zeros(1,m);
if k==1
    for j=1:m
        if x(i)<=t(j) && t(j)<x(i+1)
            y(j)=1;
        else
            y(j)=0;
        end    
    end
    return; 
end    

if abs(x(i+k)-x(i))>eps
    alpha=(t-x(i))./(x(i+k)-x(i));
else
    error('x(i+k)=x(i)')
end

if abs(x(i+k+1)-x(i+1))>eps
    beta=(x(i+k+1)-t)./(x(i+k+1)-x(i+1));
else
    error('x(i+k+1)=x(i)')
end

% vector form recursion, which, however, causes something 
% undesirable: each elements does the recursion!
y=alpha.*BSbaseVec(i,k-1,x,t)+beta.*BSbaseVec(i+1,k-1,x,t);

end

