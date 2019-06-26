% Householder Matrix : H=I-2vv'. (v is a unit vector)
% H*x=y, where x is a column vector and satisfies:||x||=||y||, y=[-alpha1,0, ... ,0].'
% Or x*H=y, where x is a row vector and satisfies:||x||=||y||, y=[-alpha1,0, ... ,0].
% The Householder matrix has the following properties:
% #1. it is Hermitian: H=H',
% #2. it is unitary: inv(H)=H',
% #3. hence it is involutory: H^2=I,
% #4. 1 is an eigenvalue of multiplicity n-1, and -1 is an eigenvalue with multiplicity 1.
% #5. det(H)=1^(n-1)*(-1)=-1.
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

alpha1=s*norm(x,2);
x(1)=x(1)+alpha1;
normx=norm(x,2);
if abs(normx)>eps
    if n==1         % column vector
        H=eye(m)-2/normx^2*x*x';
    else            % row vector
        H=eye(n)-2/normx^2*x'*x;
    end
else 
    error('Error,norm(x,2)<=eps!')
end

end
