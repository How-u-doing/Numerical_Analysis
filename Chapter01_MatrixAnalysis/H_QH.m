% QH decomposition by using Householder transformation(reflection) 
% 10170437 Mark Taylor

function [Q,H]=H_QH(A)
% H is a Hessenberg matrix

[m,n] = size(A);
if m ~= n 
    error('Input must be a square matrix!')
end

Q=eye(n);
H=A;
for j=1:n-2  
    T=eye(n);   % Property: T'==T, T'==inv(T)                                
    T(j+1:n,j+1:n)=Householder(H(j+1:n,j));  
    
    Q=Q*T;      % Q=T(1)T(2)...T(n-2).
    H=T*H*T;    % H=T(n-2)...T(2)T(1)*A*T(1)T(2)...T(n-2),
                % In particular, H is tridiagonal when A is symmetric.   
end

end

