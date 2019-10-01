% B-spline base function
function y = BSbase(i,k,x,t)
% i,k: the B-spline base of kth order(degree of k-1) whose support lies in interval (x_i,x_(i+k))
%  x : sample points, [x_1,x_2,...,x_n]
%  t : points to be evaluated 
if k==0 
    if x(i)<=t && t<x(i+1)
        y=1;
    else
        y=0;
    end
    return;
end

if x(i+k)-x(i)==0 
    alpha=0;
else
    alpha=(t-x(i))/(x(i+k)-x(i));
end

if x(i+k+1)-x(i+1)==0
     beta=0;
else
     beta=(x(i+k+1)-t)/(x(i+k+1)-x(i+1));
end

y=alpha*BSbase(i,k-1,x,t)+beta*BSbase(i+1,k-1,x,t);
end

