% Interpolation for PQRR
function u=PQRR_intpol(uc,xq,a,b)
N=length(uc)/2;   % No. of coefs of integer points
m=length(xq);     % No. of query points
u=zeros(m,1);     % approximations at these m query points
xx=linspace(a,b,N+1);
h=(b-a)/N;

phi1 = @(k)(2*k-1).*(k-1);
phi2 = @(k)4*k.*(1-k);
phi3 = @(k)k.*(2*k-1);

% a lazzy solution
for i=1:N
    for j=1:m
        if xq(j)>=xx(i) && xq(j)<=xx(i+1)
            if i==1
                u(j)=uc(1)*phi2((xq(j)-xx(1))/h)+...
                        uc(2)*phi3((xq(j)-xx(1))/h);
            else
                u(j)=uc(2*i-1)*phi2((xq(j)-xx(i))/h)+...
                        uc(2*i)*phi3((xq(j)-xx(i))/h)+...
                        uc(2*i-2)*phi1((xq(j)-xx(i))/h);
            end
        end
    end
end

% When xq is an *ordered* further division based on xx,
% for example, xq = linspace(a,b,10*N+1), then an eager
% version can be developed like this:
%{
k = (m-1)/N;
for i=1:N
    k1 = 1+k*(i-1);
    k2 = 1+k*i;
    for j=k1:k2
        if xq(j)>=xx(i) && xq(j)<=xx(i+1)
            if i==1
                u(j)=uc(1)*phi2((xq(j)-xx(1))/h)+...
                        uc(2)*phi3((xq(j)-xx(1))/h);
            else
                u(j)=uc(2*i-1)*phi2((xq(j)-xx(i))/h)+...
                        uc(2*i)*phi3((xq(j)-xx(i))/h)+...
                        uc(2*i-2)*phi1((xq(j)-xx(i))/h);
            end
        end
    end
end
%}
end

