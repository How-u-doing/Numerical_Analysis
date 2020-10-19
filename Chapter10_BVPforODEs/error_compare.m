% Error compare between PLRR and PQRR
% See error analysis at <error_analysis.mlx>
% p.w. linear: O(h^2);   cubic spline: O(h^4)
% So we can guess quadratic element is O(h^3)

% ===============================================================
% Solve a linear two-point boundary-value problem
%                   -(p(x)u'(x))'+q(x)u(x)=f(x),    a<=x<=b,
% with initial conditions
%                           u(a)=0, u'(b)=0.
% ===============================================================

% ***************************************************************
% Test 1
f=@(x)pi^2/4*sin(pi/2*x);
u=@(x)sin(pi/2*x);

% Test 2
% f=@(x)pi^2*sin(pi*x);
% u=@(x)sin(pi*x)+pi*x;
% ***************************************************************

p=@(x)x.^0; % p(x)=1
q=@(x)x.*0; % q(x)=0
a=0;
b=3;


% error = Rh = O(h^n), log(Rh)=n*log(h), slope k=n...
% which represents its convergence order. 
% So we need to let log(h) be the x axis, log(err) be the... 
% y axis and see how the slopes of these methods look like.
n=5;
N=8;    % no. of points of log(h)
h=zeros(N,1);
err=zeros(N,2);
for i=1:N  
    h(i)=(b-a)/n;
    tt=linspace(a,b,5*n+1).'; % interpolation every 5 points
    
    uc1=PLRR(f,p,q,a,b,n);
    y1=PLRR_intpol(uc1,tt,a,b);
    err(i,1)=max(abs(u(tt)-y1));
    
    uc2=PQRR(f,p,q,a,b,n);
    y2=PQRR_intpol(uc2,tt,a,b);
    err(i,2)=max(abs(u(tt)-y2));
    
    n=2*n; % delta(log(h)) = -log(2) = -0.6931
end

logh=log(h);
loge=log(err);
figure
plot(logh,loge(:,1),'b--',logh,loge(:,2),'m--','LineWidth',2)
hold on
plot(logh,2*logh,logh,3*logh)
legend('PLRR','PQRR','2*log(h)','3*log(h)','Location','northwest')
xlabel('log(h)')
ylabel('log(err)')
title('Error Comparison')

% Get convergence orders of these methods (i.e. their slopes)
% Or use built-in func polyfit instead
P1=lsq(logh,loge(:,1),1);     % use linear polynomial to fit, i.e. y = a1*x + a0, and P1 = [a0, a1] 
k_PLRR = P1(2)                % its convergence order should amount to its slope
                              % it may not be an integer, but it is pretty
                              % much close to one, say 1.9961 or whatever

P2=lsq(logh,loge(:,2),1);
k_PQRR = P2(2)


