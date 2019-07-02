% Solve nonlinear systems of equations of the form f(x)=0 via Newton's
% mothod, converting the form f(x)=0 into x=phi(x)=x-f(x)/f'(x).
% 10170437 Mark Taylor

function [x,k]=Newton(fun,x0,tol,maxit)

% Set default input arguments
if nargin<4
    maxit=1000;
    if nargin<3
        tol=1e-6;
        if nargin<2
            error('Too few input arguments!')
        end
    end
end

% Run this funtion may not as fast as bisection method, since we need to 
% call "syms" family and "matlabFunction", which demand a bit of time.
syms phi(x)
phi(x)= x-fun(x)/diff(fun,x);
phi=matlabFunction(phi);
[x,k]=fixedpoint(phi,x0,tol,maxit);

end

