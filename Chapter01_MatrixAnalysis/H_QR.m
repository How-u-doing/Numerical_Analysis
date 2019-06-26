% Householder QR Factorization
% 10170437 Mark Taylor

function [Q,R]=H_QR(A)
[m,n] = size(A);
if m ~= n 
    error('Input must be a square matrix!')
end
R=A;
Q=eye(n);
for j=1:n-1                             %     _             _
    H=eye(n);                           % H= |  Ij-1     0   | ¡ú(j-1) rows
    H(j:n,j:n)=Householder(R(j:n,j));   %    |_  0      Hj  _| ¡ú(n-(j-1)) rows
    
    R=H*R;                              % R=H(n-1)...H(2)H(1)*A                                       
    Q=Q*H;                              % Q=H(1)H(2)...H(n-1)
end

end
