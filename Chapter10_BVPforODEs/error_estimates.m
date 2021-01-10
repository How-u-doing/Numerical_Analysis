% Error Estimates for PLRR and PQRR
% See error analysis at <error_analysis.mlx>

% (1)  infinity norm --> infinity norm
% ||u-u_h||_{L^{inf}}=sqrt(int((u-u_h)^2,a,b))
% p.w. linear: O(h^2);   cubic spline: O(h^4)
% So we can guess quadratic element is O(h^3)

% (2)  L^2 norm --> L^2 norm
% ||u-u_h||_{L^2}=sqrt(int((u-u_h)^2,a,b))
% p.w. linear: O(h^2);   cubic spline: O(h^4)
% So we can guess quadratic element is O(h^3)

% (3)  H^1 norm --> L^2 norm
% ||u-u_h||_{H^1}=sqrt(int((u'-u_h')^2,a,b)+int((u-u_h)^2,a,b))
% p.w. linear: O(h);     cubic spline: O(h^3)
% So we can guess quadratic element is O(h^2)


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
u1=@(x)pi/2*cos(pi/2*x); % direvative of u

% Test 2
% f=@(x)pi^2*sin(pi*x);
% u=@(x)sin(pi*x)+pi*x;
% u1=@(x)pi*(1+cos(pi*x)); % direvative of u
% ***************************************************************

p=@(x)x.^0; % p(x)=1
q=@(x)x.*0; % q(x)=0
a=0;
b=3;


% err = O(h^n), log(err)=n*log(h), slope k = n, which represents
% its convergence order. So we need to let log(h) be the x axis,
% log(err) be the y axis, and see how the slopes of those norms
% look like.
n=5;
N=8;    % no. of points of log(h)
h=zeros(N,1);
err_inf=zeros(N,2);  % infinity norm
err_L2 = zeros(N,2); % L^2 norm
err_H1 = zeros(N,2); % H^1 norm
for i=1:N  
    h(i)=(b-a)/n;
    t =linspace(a,b,n+1).';
    tt=linspace(a,b,5*n+1).'; % interpolation every 5 points
        
    uc1=PLRR(f,p,q,a,b,n); % try also 2*n parts. See footnote for analysis
    uc2=PQRR(f,p,q,a,b,n);
    [y1, y1_1]=PLRR_intpol(t,[0;uc1],tt);
    [y2, y2_1]=PQRR_intpol(t,[0;uc2],tt);
    
    e = (u(tt)-y1).^2;
    e1 = (u1(tt)-y1_1).^2;    
    err_inf(i,1)=max(abs(u(tt)-y1)); % norm(u(tt)-y1, inf)
    err_L2(i,1)=sqrt(sum(e(1:end-1))*h(i)/5); % approx equal to:
                                              % norm(u(tt)-y1)*sqrt(h(i)/5)
    err_H1(i,1)=sqrt(sum(e(1:end-1)+e1(1:end-1))*h(i)/5);
    
    e = (u(tt)-y2).^2;
    e1 = (u1(tt)-y2_1).^2;
    err_inf(i,2)=max(abs(u(tt)-y2));
    err_L2(i,2)=sqrt(sum(e(1:end-1))*h(i)/5);
    err_H1(i,2)=sqrt(sum(e(1:end-1))*h(i)/5+sum(e1(1:end-1))*h(i)/5);
    
    n=2*n; % delta(log(h)) = -log(2) = -0.6931
end

logh=log(h);
loge_inf=log(err_inf);
loge_L2=log(err_L2);
loge_H1=log(err_H1);

% L^{inf} norm
figure
plot(logh,loge_inf(:,1),'b--',logh,loge_inf(:,2),'m--','LineWidth',2)
hold on
plot(logh,2*logh,logh,3*logh)
legend('PLRR','PQRR','2*log(h)','3*log(h)','Location','northwest')
xlabel('log(h)')
ylabel('log(err)')
title('Error Estimates under $L^\infty$ Norm','interpreter','latex')

% Get convergence orders of these methods (i.e. their slopes)
% Or use built-in func polyfit instead
P1_inf=lsq(logh,loge_inf(:,1),1);     % use linear polynomial to fit, i.e. y = a1*x + a0, and P1 = [a0, a1] 
k_PLRR_inf = P1_inf(2)                % its convergence order should amount to its slope
                                      % it may not be an integer, but it is pretty
                                      % much close to one, say 1.9961 or whatever

P2_inf=lsq(logh,loge_inf(:,2),1);
k_PQRR_inf = P2_inf(2)

% L^2 norm 
figure
plot(logh,loge_L2(:,1),'b--',logh,loge_L2(:,2),'m--','LineWidth',2)
hold on
plot(logh,2*logh,logh,3*logh)
legend('PLRR','PQRR','2*log(h)','3*log(h)','Location','northwest')
xlabel('log(h)')
ylabel('log(err)')
title('Error Estimates under $L^2$ Norm','interpreter','latex')

P1_L2=lsq(logh,loge_L2(:,1),1);
k_PLRR_L2 = P1_L2(2)

P2_L2=lsq(logh,loge_L2(:,2),1);
k_PQRR_L2 = P2_L2(2)

% H^1 norm 
figure
plot(logh,loge_H1(:,1),'b--',logh,loge_H1(:,2),'m--','LineWidth',2)
hold on
plot(logh,logh,logh,2*logh)
legend('PLRR','PQRR','log(h)','2*log(h)','Location','northwest')
xlabel('log(h)')
ylabel('log(err)')
title('Error Estimates under $H^1$ Norm','interpreter','latex')

P1_H1=lsq(logh,loge_H1(:,1),1);
k_PLRR_H1 = P1_H1(2)

P2_H1=lsq(logh,loge_H1(:,2),1);
k_PQRR_H1 = P2_H1(2)

%%%%%%%%%%%%% footnote %%%%%%%%%%%%%
% err = O(h^2) = C*h^2
% log(err') = log(C*(h/2)^2)
%           = log(err) - log(4)
% So doubling the no. of intervals
% just makes 'PLRR' translate down,
% remaining its slope unchanged.
%%%%%%%%%%%%% footnote %%%%%%%%%%%%%
