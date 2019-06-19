% Test for 
% 10170437 Mark Taylor

A=diag(15*ones(5,1))+diag(1*ones(4,1),1)+diag(1*ones(4,1),-1) +...
diag(7*ones(3,1),2)+diag(7*ones(3,1),-2) +diag(4*ones(2,1),3)+...
diag(4*ones(2,1),-3)

eg=eig(A)% to see all eigenvalues

% To see how the eigenvalue that should be closest to 17 varies as tol changes
[u9,x9,k9]=Invpower(A,17,ones(5,1),1e-9,100)
[u10,x10,k10]=Invpower(A,17,ones(5,1),1e-10,100)

% To see how the eigenvalue that should be closest to 8 varies as tol changes
[u7,x7,k7]=Invpower(A,8,ones(5,1),1e-7,100)
[u8,x8,k8]=Invpower(A,8,ones(5,1),1e-8,100)

m=(eg(4)+eg(5))/2
[u,x,k]=Invpower(A,m+.02,ones(5,1),1e-6,1000)
[u,x,k]=Invpower(A,m-.02,ones(5,1),1e-6,1000)
[u,x,k]=Invpower(A,m+.1,ones(5,1),1e-6,1000)
[u,x,k]=Invpower(A,m-.1,ones(5,1),1e-6,1000)
% Notice the changes of iterative times caused by a samll variation of q
[u,x,k]=Invpower(A,m+.2,ones(5,1),1e-6,1000)
[u,x,k]=Invpower(A,m+.2,ones(5,1),1e-6,1000)


