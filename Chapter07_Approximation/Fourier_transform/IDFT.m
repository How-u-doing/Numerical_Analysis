% Inverse Discrete Fourier Transform
function A=IDFT(x)
% x: input data vector which usually is complex type
% A: IDFT of x

% insure x is column vector
if isrow(x)
    x=x.';
end
N=length(x);
k=0:N-1;
n=0:N-1;
w=exp(1i*2*pi/N);
nk=n.'*k;
wnk=w.^(nk);
A=wnk*x/N;

end

