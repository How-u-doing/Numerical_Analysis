% Use Givens & Householder Methods to solve equation A*x=b,and return those two solutions
% 10170437 Mark Taylor

function [G_x,H_x]=G_H_x(A,b)
    % A=QR, Ax=b -->Q.'*A=R, Q.'*A*x=R*x=Q.'*b -->x=R\Q.'*b
    
    [G_Q,G_R]=G_QR(A);
    [H_Q,H_R]=H_QR(A);
    G_b=G_Q.'*b;
    H_b=H_Q.'*b;
    n=length(b);
    
    % solve R*x=b, where R is an upper triangular matrix
    G_x=zeros(n,1); G_x(n)=G_b(n)/G_R(n,n);
    H_x=zeros(n,1); H_x(n)=H_b(n)/H_R(n,n);
    for i=n-1:-1:1
        G_x(i)=(G_b(i)-G_R(i,i+1:n)*G_x(i+1:n))/G_R(i,i); 
        H_x(i)=(H_b(i)-H_R(i,i+1:n)*H_x(i+1:n))/H_R(i,i);        
    end
    
end
