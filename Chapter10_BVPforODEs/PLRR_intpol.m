% Interpolation for PLRR
function u=PLRR_intpol(uc,xq,a,b)
N=length(uc);   % No. of coefs of approx solution
m=length(xq);   % No. of query points
u=zeros(m,1);   % approximations at these m query points
xx=linspace(a,b,N+1);
h=(b-a)/N;

% a lazzy solution
for i=1:N
    for j=1:m
        if xq(j)>=xx(i) && xq(j)<=xx(i+1)
            if i==1
                u(j)=uc(1)*(xq(j)-xx(1))/h;
                continue
            end
            u(j)=uc(i-1)*(-(xq(j)-xx(i+1))/h)+uc(i)*(xq(j)-xx(i))/h;
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
                u(j)=uc(1)*(xq(j)-xx(1))/h;
                continue
            end
            u(j)=uc(i-1)*(-(xq(j)-xx(i+1))/h)+uc(i)*(xq(j)-xx(i))/h;
        end
    end
end
%}
end

