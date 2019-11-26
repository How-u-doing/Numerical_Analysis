% Romberg Simpson Integration (iterative method)
function [S,k]=Romberg_Simpson(f,a,b,maxIt,tol)
% set default args
if nargin<5
    tol=1e-6;
    if nargin<4
        maxIt=30;   % 2^30 intervals at most
        if nargin<3
            error('Too few arguments')
        end
    end
end

% R1, R2 are used to store the most recent rows
R1=zeros(1,maxIt+1);
R1(1)=CompositeInt(f,a,b,1,"Simpson");
for k=1:maxIt
    
    R2=zeros(1,maxIt+1);
    R2(1)=CompositeInt(f,a,b,2^k,"Simpson");
    for j=2:k+1
        R2(j)=1/(4^(j-1)-1)*(4^(j-1)*R2(j-1)-R1(j-1));
    end
    
    % tell when to stop
    if abs(R2(k+1)-R1(k))<tol
        S=R2(k+1);
        return
    end
    
    % update R1
    R1=R2;
end
% Maximum number of iterations exceeded (the procedure was successful).
fprintf('Maximum no. of iterations exceeded, the results of last iteration are as follows:\n')
S=R2(maxIt+1);
end

