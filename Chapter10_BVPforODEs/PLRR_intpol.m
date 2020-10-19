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

end

