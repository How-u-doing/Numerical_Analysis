% Bernstein Base function approximation test
% 10170437 Mark Taylor
f=@(x) 100.*x.*(x-0.2).*(x-0.5).*(x-0.7).*(x-0.9);
a=0;b=1;
x=a:.01:b;
xx=a:.01:b;
fxx=f(xx);
for n=0:5:50
    P=BernsteinBase(f,n,x);
    % I couldn't find an appropriate way to keep f(x) & update (Bnf)(x)
    plot(xx,fxx,'g',x,P,'--r');
    legend('Original Function','Bernstein Polynomial')
    title('Bernstein Approx.');
    grid on
    pause(0.5);
end

