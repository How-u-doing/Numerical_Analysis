% ODE test of several methods
% Another example
clear; clc; close all;
f=@(x,y)-2*x*y^2;
ya=1;   % Solution: y=1/(1+x^2);
a=0;
b=1;
n=10;
[t, y1]=Forward_Euler(f,a,b,ya,n);
[~, y2]=Backward_Euler(f,a,b,ya,n);
[~, y3]=Modified_Euler(f,a,b,ya,n);
[~, y4]=Midpoint(f,a,b,ya,n);
[~, y5]=Runge_Kutta(f,a,b,ya,n);

y=1./(1+t.^2);  % Solution

plot(t,y,'g',t,y1,'b--',t,y2,'r--',t,y3,'c--',t,y4,'m--',t,y5,'k--')
legend('Original Function','Forward Euler','Backward Euler','Modified Euler','Midpoint','Runge-Kutta')

% compare Modified Euler(explicit Trapezoidal), Trapezoidal (implicit), Midpoint(Runge-Kutta method of order 2)
% This problem, however, differs from Example_2_1, we have fitting
% effect: Modified Euler > Trapezoidal for most parts, we can see
% that implicit rule is not necessarily better than explicit, though
% more calculation was required by implicit one.
figure(2)
[t, y6]=Trapezoidal(f,a,b,ya,n);
plot(t,y,'g',t,y3,'b--',t,y4,'r--',t,y6,'c--')
legend('Original Function','Modified Euler','Midpoint','Trapezoidal')


% error = Rh = O(h^n), log(Rh)=n*log(h), slope k=n...
% which represents its convergence order. 
% So we need to let log(h) be the x axis, log(err) be the... 
% y axis and see how the slopes of these methods look like.
n=10;
N=8;    % no. of points of log(h)
h=zeros(1,N);
err=zeros(5,N);
for j=1:N
    Y=zeros(5,n+1);
    [t,Y(1,:)]=Forward_Euler(f,a,b,ya,n);
    [~,Y(2,:)]=Backward_Euler(f,a,b,ya,n);
    [~,Y(3,:)]=Modified_Euler(f,a,b,ya,n);
    [~,Y(4,:)]=Midpoint(f,a,b,ya,n);
    [~,Y(5,:)]=Runge_Kutta(f,a,b,ya,n);    
    
    h(j)=(b-a)/n;
    y=1./(1+t.^2);

    for i=1:5
        err(i,j)=max(abs(y-Y(i,:)));
    end
    
    n=2*n; % delta(log(h)) = -log(2) = -0.6931
end

logh=log(h);
loge=log(err);
figure(3)
plot(logh,loge(1,:),'b--',logh,loge(2,:),'r--',logh,loge(3,:),'c--',logh,loge(4,:),'k--',logh,loge(5,:),'m--','LineWidth',2)
hold on
plot(logh,logh,logh,2*logh,logh,4*logh)
legend('Forward Euler','Backward Euler','Modified Euler','Midpoint','Runge-Kutta','log(h)','2*log(h)','4*log(h)')
xlabel('log(h)')
ylabel('log(err)')
title('Error Comparison')

% get convergence orders of these methods (i.e. their slopes)
P1=lsq(logh,loge(1,:),1);     % use linear polynomial to fit, i.e. y = a1*x + a0, and P1 = [a0, a1] 
p_Forward_Euler = P1(2)       % its convergence order should amount to its slope
                              % it may not be an integer, but it is pretty
                              % much close to one, say 2.0177 or whatever

P2=lsq(logh,loge(2,:),1);
p_Backward_Euler = P2(2)

P3=lsq(logh,loge(3,:),1);
p_Modified_Euler = P3(2)

P4=lsq(logh,loge(4,:),1);
p_Midpoint = P4(2)

P5=lsq(logh,loge(5,:),1);
p_Runge_Kutta = P5(2)

