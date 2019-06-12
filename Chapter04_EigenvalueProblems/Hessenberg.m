% Construct upper-Hessenberg matrix via Householder transformation(reflection) 
% 10170437 Mark Taylor

function [H,Q]=Hessenberg(A)
% H is an upper-Hessenberg matrix, satisfying H=Q'*A*Q.

[m,n] = size(A);
if m ~= n
    error('Input must be a square matrix!')
end

H=A;
Q=eye(n);
% Return [H,Q] as fast as possible if A is upper-Hessenberg.
% One simple approach available in MATLAB is:
% isequal(tril(A(3:n,1:n-2)),zeros(n-2))==true, which requires (n-2)^2 
% flops of compare and maybe an extra memory consumption(zeros(n-2)).
% We adopt follows steps in order to convert it into a c++ source file.

% /*****
isZero=true;    % flag that tells if A is upper-Hessenberg 
for i=3:n
    for j=1:i-2
        if abs(A(i,j))>eps
            isZero=false;
            break;
        end
    end
    if isZero==false
        break;
    end
end

if isZero==true
    return;
end
% *****/

for j=1:n-2  
    T=eye(n);   % Property: T'==T, T'==inv(T)                                
    T(j+1:n,j+1:n)=Householder(H(j+1:n,j));  
    
    Q=Q*T;      % Q=T(1)T(2)...T(n-2).
    H=T*H*T;    % H=T(n-2)...T(2)T(1)*A*T(1)T(2)...T(n-2),
                % In particular, H is tridiagonal when A is symmetric.   
end

end

