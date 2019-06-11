% Compute eigenvalues via QR algorithm
% 10170437 Mark Taylor

function [eigval,k]=Eig_test(A,tol,N)

[m,n] = size(A);
if m ~= n 
    error('A must be square!')
end

% Set default function input arguments
if nargin<3
    N=100;
    if nargin<2
        tol=1e-6;
        if nargin<1
            error('Too few input argument(s)')
        end
    end
    
end


% Construct upper-Hessenberg matrix
H=Hessenberg(A); % H is a Hessenberg matrix who has identical eigenvalues of A

eigval=zeros(n,1);
t=1;
k=1;
while k<=N
    if abs(H(n,n-1))<tol
        eigval(t)=H(n,n);t=t+1;        
        n=n-1;
    end
    
    if abs(H(2,1))<tol
        eigval(t)=H(1,1);t=t+1;
        n=n-1;
        H(1,1)=H(2,2);
        for j=2:n
            H(j,j)=H(j+1,j+1);
            H(j-1,j)=H(j,j+1);
            H(j,j-1)=H(j+1,j);
        end
    end
        
    if n==1
        eigval(t)=H(1,1);
        return;
    end
       
% /***** Compute shift
    b=-(H(n-1,n-1)+H(n,n));
    c=H(n-1,n-1)*H(n,n)-H(n-1,n)*H(n,n-1);
    d=sqrt(b^2-4*c);
    u1=(d-b)/2;
    u2=-(b+d)/2;
    
    if n==2
        eigval(t)=u1;t=t+1;
        eigval(t)=u2;
        return;
    end
    
    % Choose shift that is closer to H(n,n)
    if abs(u1-H(n,n))<abs(u2-H(n,n))
        s=u1;
    else
        s=u2;
    end
     
%  *****/ 


    % Perform shift
    % H=H-s*I
    for i=1:n
        H(i,i)=H(i,i)-s;
    end
    
    % Compute H_(k+1) via Givens rotation
    for j=1:n-1
        P = GV(H(j,j),H(j+1,j),j,j+1,n);
        % H=RQ=P(n-1)...p(2)P(1)*H*P(1)'P(2)'...P(n-1)'       
        H(1:n,1:n)=P*H(1:n,1:n)*P.';        
        %/*
        % Indeed, row 83 and 85 can be optimized to cut down the amount of 
        % calculation due to the fact that P left(right) multiplies H only 
        % changes two rows(columns) of H. See Eig.m for implementation.
        %*/
    end
    % H=H+s*I
    for i=1:n
        H(i,i)=H(i,i)+s;
    end    
    
    k=k+1; 
end

% The number of iterations was exceeded.
k=k-1;% k=N
fprintf('\nCannot compute the approximate eigenvalues within %d iterations in the tolerance of %d!\n',N,tol);
if t>1
    fprintf('Calculated partial eigenvalues in the last iteration are as followed:\n');
    eigval=eigval(1:t-1);
else
    eigval='Didn''t figure out one yet';
end

end