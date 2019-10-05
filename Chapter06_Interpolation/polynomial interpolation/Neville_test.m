% Neville's method test
% 10170437 Mark Taylor
%% test on f(x)= 1/(1+x^2)
figure('Name','1/(1+x^2) on [-5,5]');
x=-5:5;
y=1./(1+x.^2);
xq=-5:.1:5;
s = Neville(x,y,xq);
s2=NI(x,y,xq);
% Indeed, Newton's & Neville's interpolating polynomials are identical
% verify the fact from the figure
plot(x,y,'o',xq,s,'.r',xq,s2,'--b',xq,1./(1+xq.^2),'g');
legend('Sample Points','Neville''s','Newton''s','Original Function')
title('Neville''s Interpolation')
grid on

figure('Name','1/(1+x^2) on [-7,7]');
x=-7:7;
y=1./(1+x.^2);
xq=-7:.1:7;
s = Neville(x,y,xq);
plot(x,y,'o',xq,s,'--r',xq,1./(1+xq.^2),'g');
legend('Sample Points','Neville''s','Original Function')
title('Neville''s Interpolation')
grid on

figure('Name','1/(1+x^2) on [-10,10]');
x=-10:10;
y=1./(1+x.^2);
xq=-10:.1:10;
s = Neville(x,y,xq);
plot(x,y,'o',xq,s,'--r',xq,1./(1+xq.^2),'g');
legend('Sample Points','Neville''s','Original Function')
title('Neville''s Interpolation')
grid on


%% test on f(x)= sin(x)
figure('Name','Perfect approximation on sin(x)');
x=-10:10;
y=sin(x);
xq=-10:.1:10;
s = Neville(x,y,xq);
plot(x,y,'o',xq,s,'.r',xq,sin(xq),'g');
legend('Sample Points','Neville''s','Original Function')
title('Neville''s Interpolation')
grid on

%% test on f(x)= x^3-x
figure('Name','x^3-x');
x=-10:10;
y=x.^3-x;
xq=-10:.1:10;
s = NI(x,y,xq);plot(x,y,'o',xq,s,'.r',xq,xq.^3-xq,'g');
legend('Sample Points','Neville''s','Original Function')
title('Neville''s Interpolation')
grid on
