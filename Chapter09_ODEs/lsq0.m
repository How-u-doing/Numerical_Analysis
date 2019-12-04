% Least Squares Approximation via Span{1,x,x^2,...,x^n}
function a=lsq0(x,y,n,rho)
%INPUT:
    % (x,y): pairs of data, vectors
    % n: approximate polynomial of degree at most n
    % rho: weights, vector
%OUPUT:
    % a: coefficients in row vector form
    
% insure x, y & rho are column vector
if isrow(x)
    x=x.';
end
if isrow(y)
    y=y.';
end
N=length(x);    % no. of samples
C=zeros(N,n+1);
C(:,1)=ones(N,1);
for j=2:n+1
    C(:,j)=x.*C(:,j-1);
end
% fast return if weights are even
if nargin<4
    a=(C.'*C)\(C.'*y);
    a=a.';
    return;
elseif isrow(rho)
    rho=rho.';
end

A=zeros(n+1);
Y=zeros(n+1,1);
for i=1:n+1
    for j=1:n+1
        A(i,j)=sum(C(:,i).*C(:,j).*rho);        
    end
    Y(i)=sum(C(:,i).*y.*rho);
end
% Solve normal equation: A*X=Y, X=[a0 a1 ... an].'
a=A\Y;
a=a.';
end

