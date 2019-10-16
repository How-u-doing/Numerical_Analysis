% To find a solution to f(x) = 0 given the continuous function f on
% the interval [x0,x1], where f(x0) and f(x1) have opposite signs.
% Note that this method is NOT as fast as scant method, however,
% this method can control all iterative zeros within interval [x0,x1],
% whereas the latter doesn't guarantee it.
% 10170437 Mark Taylor
function [x,k]= FalsePosition(f,x0,x1,tol,maxIt)
% Set default input arguments
if nargin<5
    maxIt=1000;
    if nargin<4
        tol=1e-6;
        if nargin<3
            error('Too few input arguments!')
        end
    end
end
fx0=f(x0);fx1=f(x1);
if sign(fx0)==sign(fx1)
    error('f(x0)*f(x1)>=0! Please insure f(x0)*f(x1)<0 to guarantee there exists at least one zero over£¨a,b)!')
end
for  k=1:maxIt
    x=x0-(x0-x1)/(fx0-fx1)*fx0;
    fx=f(x);
    if abs(fx)<tol
        return;
    end
    if sign(fx0)~=sign(fx)
        x1=x;
        fx1=fx;
    else
        x0=x;
        fx0=fx;
    end
end
% The number of iterations was exceeded.
fprintf('\nCannot compute the approximate solution within %d iterations in the tolerance of %d!\n',maxIt,tol);
fprintf('Calculated approximate solution in the last iteration is as followed:\n');
end

