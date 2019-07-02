% Solve nonlinear systems of equations of the form f(x)=0 via fixed point mothod.
% 10170437 Mark Taylor

function [x,k]=fixpoint(phi,x0,tol,maxit)
% x=phi(x), 


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


k=0;                                            % Practical iterative times
x=x0;
while k<maxit
    if abs(x-phi(x))<tol
        return;
    end
    
    x=phi(x);
    
    k=k+1;
end
% The number of iterations was exceeded.
k=k-1;                                          % k=maxit
fprintf('\nCannot compute the approximate solution within %d iterations in the tolerance of %d!\n',maxit,tol);
fprintf('Calculated solution in the last iteration is as followed:\n');
end

