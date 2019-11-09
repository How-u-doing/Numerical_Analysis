% Inverse Fast Fourier Transform
function A=My_IFFT(x)
% x: input data vector, which usually is complex type
% A: N  Inverse Discrete Fourier Transform (IDFT) points
N=length(x);    % N=2^p, N usually is very large
p=log2(N);
% tell if log2(N) is an integer to proceed
if mod(p,2)~=0 && mod(p,2)~=1
    error('log2(N) must be an integer!')
end
w=0:N/2-1;
w=exp(1i*2*pi/N).^w;
if isrow(x)
    x1=zeros(1,N);
else
    x1=zeros(N,1);
end
for q=1:p
    if(mod(q,2)==1)
        for k=0:2^(p-q)-1
            for j=0:2^(q-1)-1
                x1(k*2^q+j+1)=x(k*2^(q-1)+j+1)+x(k*2^(q-1)+j+2^(p-1)+1);
                x1(k*2^q+j+2^(q-1)+1)=(x(k*2^(q-1)+j+1)-...
                    x(k*2^(q-1)+j+2^(p-1)+1))*w(k*2^(q-1)+1);
            end            
        end
    else
        for k=0:2^(p-q)-1
            for j=0:2^(q-1)-1
                x(k*2^q+j+1)=x1(k*2^(q-1)+j+1)+x1(k*2^(q-1)+j+2^(p-1)+1);
                x(k*2^q+j+2^(q-1)+1)=(x1(k*2^(q-1)+j+1)-...
                    x1(k*2^(q-1)+j+2^(p-1)+1))*w(k*2^(q-1)+1);
            end
        end
    end
end
if mod(p,2)==0
    % p is even
    A=x/N;
else % p is odd
    A=x1/N;
end

end

