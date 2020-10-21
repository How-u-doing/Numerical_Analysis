% Interpolation for PLRR
function u=PLRR_intpol(uc,xq,a,b)
N=length(uc);   % No. of coefs of approx solution
m=length(xq);   % No. of query points
u=zeros(m,1);   % approximations at these m query points
xx=linspace(a,b,N+1);
h=(b-a)/N;

% a lazzy solution
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

% When xq is an *ordered* further division based on xx,
% for example, xq = linspace(a,b,10*N+1), then an eager
% version can be developed like this:
% {
k = (m-1)/N;
for i=1:N
    % Note that, we are better off using floor, since floor(k * N)
    % might NOT be equal to m-1, in turn, xq(m) remains what it is.
    k1 = 1 + ceil(k*(i-1)); % Get upper bounds, at most 1 element of
    k2 = 1 + ceil(k*i);     % xq(k1:k2) NOT in [xx(i), xx(i+1)].
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

% ************************FOOTNOTE***************************
% Say, m=101, N=3, floor(k*N)=100?? Luckily, this is true in 
% MATLAB (C, C++, Python, etc.) due to the fact that: 
% >> k=100/3; % k = 33.333333333333336
% >> k*3==100
% ans = (logical) 1
% That's how double-precision numbers store and operate.
% Notice that floor(33.3333*3) = 99, NOT 100.
% ************************FOOTNOTE***************************
end

