% Least Squares Approximation via Span{1,x,x^2,...,x^n}
function a=lsq(x,y,n,rho)
%INPUT:
    % (x,y): pairs of data, vectors
    % n: approximate polynomial of degree at most n
    % rho: weights, vector
%OUPUT:
    % a: coefficients of row vector form in ASCENDING order
    
% insure x, y & rho are column vector
if isrow(x)
    x=x.';
end
if isrow(y)
    y=y.';
end

N=length(x);    % no. of samples
if nargin<4
    rho=ones(N,1);
elseif isrow(rho)
    rho=rho.';
end

A=zeros(N,n+1); % sqrt weighted coefficient matrix
sr=sqrt(rho);   % when N is enormous£¬it may consume some time,
                % which, however, can be eluded. (see lsq0.m)
A(:,1)=ones(N,1).*sr;
for j=2:n+1
    % A(:,j)=sr.*x.^(j-1); % may be slow (when n is very large) 
    A(:,j)=x.*A(:,j-1);
end
y=sr.*y; % sqrt weighted data value
% Solve normal equation: A.'*A*X=A.'*Y, X=[a0 a1 ... an].'
a=(A.'*A)\(A.'*y);
a=a.';
end

