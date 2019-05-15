 % Cramer's Rule: solution of the system of linear equations of the form Ax=b
 % 10170437 Mark Taylor
function x= Cramer(A, b)
[m,n] = size(A);
if m ~= n % make sure A is a square matrix to continue following steps
    error('Input must be a square matrix!!')
end

D=det(A);
if abs(D)<=eps 
% In freemat, '=' in '<=' is necessary, otherwise it cannot detect 
% det(A£©==0 and get resultant wrong solution of this equation.
    error('Error! det(A)==0')
end
if abs(D)<1e-4
    warning('det(A) is very tiny(less than 1e-4), great error may transpire')
end

x=zeros(n,1);
for j=1:n
    B=A;
    B(:,j)=b;
    x(j)=det(B)/D;
end

end
