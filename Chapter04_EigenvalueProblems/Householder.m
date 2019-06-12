% Householder Matrix : H=I-2vv'. (v is a unit vector)
% H*x=y, where x is a column vector and satisfies:||x||=||y|| ,y=[-¦Á1,0, ... ,0]'
% Or x*H=y,where x is a row vector and satisfies:||x||=||y|| ,y=[-¦Á1,0, ... ,0]
% 10170437 Mark Taylor

function H = Householder(x)
[m,n]=size(x);
if min(m,n)~=1
    error('Input must be a vector!')
end

if x(1)>=0
    s=1;
else
    s=-1;
end
alpha1=s*norm(x);   % alpha1=s(x)*||x||
x(1)=x(1)+alpha1;
if abs(norm(x))>eps
    if n==1         % column vector
        H=eye(length(x))-2/(norm(x))^2*x*x.';
    else            % row vector
        H=eye(length(x))-2/(norm(x))^2*x.'*x;
    end
else 
    error('Error,norm(x,2)<=eps!')
end

end
