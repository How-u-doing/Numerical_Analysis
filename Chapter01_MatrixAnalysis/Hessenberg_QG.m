% Q.'*A*Q=G, where Q is a orthogonal matrix and G is a Hessenberg matrix 
% 10170437 Mark Taylor
function [Q,G]=Hessenberg_QG(A) % Implemented by using Householder transformation(reflection)
    
[m,n] = size(A);
if m ~= n % make sure A is a square matrix to continue following steps
    error('Input must be a square matrix!')
end
copy_A=A;
Q=eye(n);
for j=1:n-2                                     %     _           _
    H=eye(n);                                   % H= |   Ij    0   | ¡új rows
    H(j+1:n,j+1:n)=Householder(A(j+1:n,j));     %    |_  0     Hj _| ¡ú(n-j) rows   
    A=H*A*H; %  A=H.'*A*H (H is symmetric)
    Q=Q*H; %  Q=H(1)H(2)...H(n-2)
end
G=Q.'*copy_A*Q; % Especially when A is symmetric, G is a tridiagonal matrix
end

