% Cubic Spline with clamped boundary
% 10170437 Mark Taylor 
function s=CubicSpline(xq,x,y,dya,dyb)
% xq: query points, specified as a vector, to be evaluated
% x,y: vector x consists of n sample points, and vector y is its corresponding function values
% dya,dyb: the derivatives f'(x_0) and f'(x_N) respectively
% s: the spline function values at those query points
k=length(xq);        % the no. of query points
n=length(x);         % the no. of sample points
if n~=length(y)
    error('mismatch between x and y\n')
end
if nargin<=3
    dya=0;  dyb=0;    
end

A=2*eye(n);         % coefficient matrix
A(1,2)=1;   A(n,n-1)=1;
d=zeros(n,1);       % right-hand side vector
d(1)=6/(x(2)-x(1))*((y(2)-y(1))/(x(2)-x(1))-dya);       % 6f[x_1,x_1,x_2]
d(n)=6/(x(n)-x(n-1))*(dyb-(y(n)-y(n-1))/(x(n)-x(n-1))); % 6f[x_(n-1),x_n,x_n]
for i=2:n-1
    h0=x(i)-x(i-1);
    h1=x(i+1)-x(i);
    A(i,i+1)=h1/(h1+h0);        % lambda_i
    A(i,i-1)=1-A(i,i+1);        % mu_i    
    d(i)=6*((y(i+1)-y(i))/h1-(y(i)-y(i-1))/h0)/(h0+h1); % 2nd divided difference * 6
end
M=solveTridiag(A,d);            % solve A*M=d

% construct cubic spline s(x) and compute those query points
s=zeros(k,1);
for i=2:n
    for j=1:k
        if xq(j)<=x(i) && xq(j)>=x(i-1)            
            x0=xq(j)-x(i-1);
            x1=x(i)-xq(j);
            h=x(i)-x(i-1);
            s(j)=(M(i-1)*x1^3+M(i)*x0^3)/6/h+(y(i-1)*x1+y(i)*x0)/h-...
                (M(i-1)*x1+M(i)*x0)*h/6;
        end
    end    
end

end

