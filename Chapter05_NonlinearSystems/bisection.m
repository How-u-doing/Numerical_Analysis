% Solve nonlinear systems of equations of the form f(x)=0 via bisection mothod.
% 10170437 Mark Taylor

function [x,k]=bisection(fun,a,b,tol,maxit)

% Set default input arguments
if nargin<5
    maxit=1000;
    if nargin<4
        tol=1e-6;
        if nargin<3
            error('Too few input arguments!')
        end
    end
end


k=0;                                            % Practical iterative times
if a>=b
    error('a must be less than b, i.e. a<b!')
elseif abs(fun(a))<tol
    x=a;
    return;
elseif abs(fun(b))<tol
    x=b;
    return;
end

r=(a+b)/2;
while k<=maxit
    if abs(fun(r))<tol
        x=r;
        return;
    end
    
    if fun(a)*fun(r)<0
        b=r;
    elseif fun(b)*fun(r)<0
        a=r;
    end
    
    r=(a+b)/2;    
    k=k+1;
end
% The number of iterations was exceeded.
k=k-1;                                          % k=maxit
x=r;
fprintf('\nCannot compute the approximate solution within %d iterations in the tolerance of %d!\n',maxit,tol);
fprintf('Calculated approximate solution in the last iteration is as followed:\n');
end

