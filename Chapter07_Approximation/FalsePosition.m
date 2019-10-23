% To find a solution to Pn(x)-f(x) = 0 given polynomial coefficents p and 
% the continuous function f on the interval [x0,x1], where Pn(x0)-f(x0) & 
% Pn(x1)-f(x1) have opposite signs or their product is zero.
% Note that this method is NOT as fast as scant method, however, it can 
% control all iterative zeros within interval [x0,x1], whereas the latter
% doesn't guarantee it.
% 10170437 Mark Taylor
function [x,k]= FalsePosition(f,p,x0,x1,tol,maxIt)
% Set default input arguments
if nargin<6
    maxIt=1000;
    if nargin<5
        tol=1e-6;
        if nargin<4
            error('Too few input arguments!')
        end
    end
end
g0=polyval(p,x0)-f(x0);
g1=polyval(p,x1)-f(x1);
if abs(g0)<=1e-6
    x=x0;k=0;
    return;
elseif abs(g1)<=1e-6
    x=x1;k=0;
    return;
elseif sign(g0)==sign(g1)
    error('Pn(x)-f(x) must change sign on endpoints')
end
for  k=1:maxIt
    x=x0-(x0-x1)/(g0-g1)*g0;
    gx=polyval(p,x)-f(x);
    if abs(gx)<tol
        return;
    end
    if sign(g0)~=sign(gx)
        x1=x;
        g1=gx;
    else
        x0=x;
        g0=gx;
    end
end
% The number of iterations was exceeded.
fprintf('\nCannot compute the approximate solution within %d iterations in the tolerance of %d!\n',maxIt,tol);
fprintf('Calculated approximate solution in the last iteration is as followed:\n');
end

