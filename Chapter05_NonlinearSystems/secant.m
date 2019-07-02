% Solve nonlinear systems of equations of the form f(x)=0 via secant mothod.
% 10170437 Mark Taylor

function [x,k]=secant(fun,x1,x2,tol,maxit)

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


k=1;                                            % Practical iterative times
while k<=maxit
    
    fx1=fun(x1);
    x=x1-(x1-x2)/(fx1-fun(x2))*fx1;
    if abs(fun(x))<tol
        return;
    end
    
    x1=x2;
    x2=x;
    
    k=k+1;
end
% The number of iterations was exceeded.
k=k-1;                                          % k=maxit
fprintf('\nCannot compute the approximate solution within %d iterations in the tolerance of %d!\n',maxit,tol);
fprintf('Calculated approximate solution in the last iteration is as followed:\n');
end

