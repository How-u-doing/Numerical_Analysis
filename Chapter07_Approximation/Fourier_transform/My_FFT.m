% Fast Fourier Transform
function x=My_FFT(A)
% A: input data vector, which can be complex type
% x: N  Discrete Fourier Transform (DFT) points
N=length(A);    % N=2^p, N usually is very large
p=log2(N);
% tell if log2(N) is an integer to proceed
if mod(p,2)~=0 && mod(p,2)~=1
    error('log2(N) must be an integer!')
end
w=0:N/2-1;
w=exp(-1i*2*pi/N).^w;
if isrow(A)
    A1=zeros(1,N);
else
    A1=zeros(N,1);
end
for q=1:p
    if(mod(q,2)==1)
        for k=0:2^(p-q)-1
            for j=0:2^(q-1)-1
                A1(k*2^q+j+1)=A(k*2^(q-1)+j+1)+A(k*2^(q-1)+j+2^(p-1)+1);
                A1(k*2^q+j+2^(q-1)+1)=(A(k*2^(q-1)+j+1)-...
                    A(k*2^(q-1)+j+2^(p-1)+1))*w(k*2^(q-1)+1);
            end            
        end
    else
        for k=0:2^(p-q)-1
            for j=0:2^(q-1)-1
                A(k*2^q+j+1)=A1(k*2^(q-1)+j+1)+A1(k*2^(q-1)+j+2^(p-1)+1);
                A(k*2^q+j+2^(q-1)+1)=(A1(k*2^(q-1)+j+1)-...
                    A1(k*2^(q-1)+j+2^(p-1)+1))*w(k*2^(q-1)+1);
            end
        end
    end
end
if mod(p,2)==0
    % p is even
    x=A;
else % p is odd
    x=A1;
end

end

