% Householder transformation &  QR factorization
% 10170437 Mark Taylor
function [Q,R]=H_QR(A) % Householder QR decomposition method
[m,n] = size(A);
if m ~= n % make sure A is a square matrix to continue following steps
    error('Input must be a square matrix!')
end
copy_A=A;
Q=eye(n);
for j=1:n-1                                     %     _             _
    H=eye(n);                                   % H= |  Ij-1     0   | ¡ú(j-1) rows
    H(j:n,j:n)=Householder(A(j:n,j));           %    |_  0      Hj  _| ¡ú(n-(j-1)) rows
    A=H*A;                                          
    Q=H*Q; %   Q=H(n-1)*H(n-2)*...*H(1)*In
end
R=Q*copy_A;
Q=Q.';
end
