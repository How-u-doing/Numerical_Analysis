% Piecewise Cubic Hermite Interpolating 
% 10170437 Mark Taylor 
function s=PCHI(xq,x,y,yd)
% xq: query points, specified as a vector, to be evaluated
% x,y: vector x & y are sample points and their corresponding function values, respectively
% yd: corresponding derivatives at those sample points
% s: the function values at those interpolating£¨query£© points
m=length(x);    % No. of sample points
n=length(xq);   % No. of query points
s=zeros(n,1);   % approximations at these n query points
for i=1:m-1
    h=x(i+1)-x(i);
    for j=1:n
        % Note that if some query point is last sample point, i.e. xq(j)==x(n) for some j,
        % then we can't calculate that one, u can change < as <=. However when query points 
        % contain all those samples points, we need calculate n-2 more duplicate points. 
        % Your call!
        if x(i)<=xq(j) && xq(j)<x(i+1)            
            x0=xq(j)-x(i);
            x1=xq(j)-x(i+1);
            s(j)=y(i)*(1+2*x0/h)*(x1/h)^2+y(i+1)*(1-2*x1/h)*(x0/h)^2+...
                yd(i)*x0*(x1/h)^2+yd(i+1)*x1*(x0/h)^2;
        end
    end   
end

end

