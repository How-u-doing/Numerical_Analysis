 % Cramer's Rule: solution of a system of linear equations Ax=b
 % 10170437 Mark Taylor
function x = Cramer(A,b)
[m,n] = size(A);
if m ~= n
    error('Input must be square!')
end
D=det(A);
if abs(D)<eps
    error('Error, det(A)=0!')
end

x=zeros(n,1);
for j=1:n
    B=A;
    B(:,j)=b;
    x(j)=det(B)/D;
end
end
