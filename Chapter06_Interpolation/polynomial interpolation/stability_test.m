% Stability test of interpolating polynomials of different degrees

t=-5:0.1:5;
y=1./(1.+t.^2);
syms x;
% construct 6 sample points
x1=(-5:2:5).';
y1=1./(1.+x1.^2);
P5=simplify(Newton([x1,y1]));
P5=vpa(P5,7);
Y1=subs(P5,x,t);

% construct 11 sample points
x2=(-5:5).';
y2=1./(1.+x2.^2);
P10=simplify(Newton([x2,y2]));
P10=vpa(P10,7);
Y2=subs(P10,x,t);

% draw an image that consists of P5(x),P10(x) & f(x)=1/(1+x^2)
figure('Name','stability test','NumberTitle','off');
plot(x2,y2,'o',t,Y1,'r',t,Y2,'b',t,y,'k');
legend('sample points','P5(x)','P10(x)','f(x)');
title('Unstability of high degree interpolating polynomial');

