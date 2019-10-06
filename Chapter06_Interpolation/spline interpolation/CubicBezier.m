% Cubic Bezier Curve
% 10170437 Mark Taylor
key=input('Press 1 for inputing control points via mouse or 2 via keyboard repectively-->');
if key==1
	axis([-100 100 -100 100])
    p=ginput(4);
else
    if key==2
        p=input('Enter 4 control points within the square brackets [x1 y1;x2 y2;x3 y3;x4 y4]:\n--> ');
    else
        error('Invalid key!')
    end
end

CubicBezierPlot(p(:,1),p(:,2))

function CubicBezierPlot(x,y)
t=0:0.01:1;
xt= (1-t).^3*x(1) + 3*(1-t).^2.*t*x(2) +...
    3*(1-t).*t.^2*x(3) + t.^3*x(4);
yt= (1-t).^3*y(1) + 3*(1-t).^2.*t*y(2) +...
    3*(1-t).*t.^2*y(3) + t.^3*y(4);
plot(x,y,'.',xt,yt,'--r',x,y,'b')
text(x,y,{'P1','P2','P3','P4'})
legend('Control Points','Bezier Curve')
end
