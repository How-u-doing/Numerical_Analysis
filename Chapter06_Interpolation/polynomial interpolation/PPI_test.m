% Piecewise Polynomial Interpolation test
% 10170437 Mark Taylor
x=-10:10;
y=1./(x.^2+1);
yd=-2.*x./(1+x.^2).^2;
xq=-10:0.1:10;
s1=PLI(xq,x,y);                         % Piecewise Linear Interpolation
s2=PCHI(xq,x,y,yd);                     % Piecewise Cubic Hermite Interpolation
s3=CubicSpline(xq,x,y,yd(1),yd(11));    % Cubic Spline
plot(x,y,'-o',xq,s1,'--r',xq,s2,'.b',xq,s3,'-.k',xq,1./(xq.^2+1),'g');
legend('Sample Points, system','PLI','PCHI','Cubic Spline','Original Function');
title('1/(1+x^2)');
% You are capable of seeing following facts from the figure(zoom in to see more clearly)
% No.1: The PLI line segments(--r) pricesely are the lines between sample points drawn by system(-o) 
% No.2: PCHI is more accurate than cubic spline, since the former require more info (all derivatives)

