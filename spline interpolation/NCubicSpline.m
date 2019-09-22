% Natural Cubic Spline (natural or free boundary)
% 10170437 Mark Taylor 
function s=NCubicSpline(xq,x,y)
% xq: query points, specified as a vector, to be evaluated
% x,y: vector x consists of n sample points, and vector y is its corresponding function values
% s: the natural spline function values at those query points
k=length(xq);        % the no. of query points
n=length(x);         % the no. of sample points
if n~=length(y)
    error('mismatch between x and y\n')
end

A=2*eye(n);         % coefficient matrix
% A(1,2)=A(n,n-1)=0
d=zeros(n,1);       % right-hand side vector
% d(1)=d(n)=0 
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

