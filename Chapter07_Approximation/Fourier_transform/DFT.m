% Discrete Fourier Transform
function x=DFT(A)
% A: input data vector which can be complex type
% x: DFT of A

% insure A is column vector
if isrow(A)
    A=A.';  % do NOT write A=A' 
end
N=length(A);
k=0:N-1;
n=0:N-1;
w=exp(-1i*2*pi/N);
nk=n.'*k;
wnk=w.^nk;  % the coef matrix
x=wnk*A;

end

