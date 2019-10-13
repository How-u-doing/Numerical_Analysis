% Bernstein Polynomial approximation test
% 10170437 Mark Taylor
f=@(x) 1./(1+x.^2);
a=-2;b=2;   % interval [a,b] has some influence on the approx effect   
x=a:.01:b;
xx=a:.01:b;
fxx=f(xx);
for n=0:5:50
    P=BernsteinPoly(f,n,x,a,b);
    % I couldn't find an appropriate way to keep f(x) & update (Bnf)(x)
    plot(xx,fxx,'g',x,P,'--r');
    legend('Original Function','Bernstein Polynomial')
    title('Bernstein Approx.');
    grid on
    pause(0.5);
end

% change the interval to see the approx effect
a=-10;b=10;
x=a:.01:b;
xx=a:.01:b;
fxx=f(xx);
figure
for n=0:5:50
    P=BernsteinPoly(f,n,x,a,b);
    plot(xx,fxx,'g',x,P,'--r');
    legend('Original Function','Bernstein Polynomial')
    title('Bernstein Approx.');
    grid on
    pause(0.5);
end