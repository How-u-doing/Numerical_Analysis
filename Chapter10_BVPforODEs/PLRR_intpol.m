% % Interpolation for PLRR
% function u=PLRR_intpol(uc,xq,a,b)
% N=length(uc);   % No. of coefs of approx solution
% m=length(xq);   % No. of query points
% u=zeros(m,1);   % approximations at these n query points
% xx=linspace(a,b,N+1);
% h=1/N;
% for i=1:N
%     for j=1:m
%         if xq(j)>=xx(i) && xq(j)<=xx(i+1)
%             if i==1
%                 u(j)=uc(1)*(xq(j)-xx(1))/h;
%                 continue
%             end
%             u(j)=uc(i-1)*(-(xq(j)-xx(i+1))/h)+uc(i)*(xq(j)-xx(i))/h;
%         end
%     end
% end
% end

% Interpolation for PLRR only at those N+1 evenly spaced points
function u=PLRR_intpol(uc,a,b)
N=length(uc);     % No. of coefs of approx solution
u=zeros(N+1,1);   % approximations at these N+1 evenly spaced points
xx=linspace(a,b,N+1);
h=1/N;
% u(0)=0;
for i=2:N
    u(i)=uc(i-1)*(-(xx(i)-xx(i+1))/h);
end
u(N+1)=uc(i)*(xx(N+1)-xx(N))/h;

end

