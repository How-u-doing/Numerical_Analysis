% Piecewise Linear Interpolation
% 10170437 Mark Taylor 
function s=PLI(xq,x,y)
% xq: query points, specified as a vector, to be evaluated
% x,y: vector x & y are sample points and their corresponding function values, respectively
% s: the function values at those interpolating£¨query£© points
m=length(x);    % No. of sample points
n=length(xq);   % No. of query points
s=zeros(n,1);   % approximations at these n query points
for i=1:m-1
    for j=1:n
        % Note that if some query point is last sample point, i.e. xq(j)==x(n) for some j,
        % then we can't calculate that one, u can change < as <=. However when query points 
        % contain all those samples points, we need calculate n-2 more duplicate points. 
        % Your call!
        if x(i)<=xq(j) && xq(j)<x(i+1)
            s(j)=(xq(j)-x(i+1))/(x(i)-x(i+1))*y(i)+(xq(j)-x(i))/(x(i+1)-x(i))*y(i+1);
        end
    end   
end

end

