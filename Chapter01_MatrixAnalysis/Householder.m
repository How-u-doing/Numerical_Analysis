% Householder Matrix : H=I-2vv'. (v is a unit vector)
% H*x=y, where x is a column vector and satisfies:||x||=||y|| ,y=[-��1,0, ... ,0]'
% Or x*H=y,where x is a row vector and satisfies:||x||=||y|| ,y=[-��1,0, ... ,0]
% 10170437 Mark Taylor
function H = Householder(x)
[m,n]=size(x);
if min(m,n)~=1 % x is neither a column vector nor a row vector ,but a matrix
    error('Input error, x must be a vector!')
end

if x(1)>=0
    s=1;
else
    s=-1;
end
alpha1=s*norm(x);% ��1=s(x)*||x||
% Here I adopt norm(x) to compute ||x|| rather than sqrt(x.'*x), whose purpose  
% is to make sure both row and column vector x can be calculated as a fixed
% number ||x||.(When x is a row vector,sqrt(x.'*x) is a n-by-n matrix,whose
% each element is sqrted)
x(1)=x(1)+alpha1;
if abs(norm(x))>eps
    if n==1    % column vector
        H=eye(length(x))-2/(norm(x))^2*x*x.';
    else       % row vector
        H=eye(length(x))-2/(norm(x))^2*x.'*x;
    end
else 
    error('Error,||x+��1e1||<=eps!')
end
end
