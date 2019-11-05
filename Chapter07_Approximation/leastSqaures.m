% Least Squares Approximation via Span{1,x,x^2,...,x^n}
function a=leastSqaures(x,y,n)
%INPUT:
    % (x,y): pairs of data, vectors
    % n: approximate polynomial of degree at most n
%OUPUT:
    % a: coefficients in row vector form
    
% insure x, y are column vector
if isrow(x)
    x=x.';
end
if isrow(y)
    y=y.';
end
N=length(x);    % no. of data points
A=zeros(N,n+1);
A(:,1)=ones(N,1);
for j=2:n+1
    % slow: A(:,j)=x.^(j-1) 
    A(:,j)=x.*A(:,j-1); % faster
end
% Solve normal equation: A.'*A*X=A.'*Y, X=[a0 a1 ... an]
% Note that as n goes large, A.'*A is usually ill-conditioned.
a=(A.'*A)\(A.'*y);
a=a.';
end

