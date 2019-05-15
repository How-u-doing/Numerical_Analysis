
% Test: perform "Aasen_test.m" without selecting pivots.
% You can load "AasenTestData.mat" which was generated from 
% "Aasen0.m"  to permutate A, that is to say Enter A=P*A*P.'  
% in your command window will generate a good condition A,
% which enables you perform elimination with no need to pivoting.  

function [P,L,T]=Aasen_test(A)

[m,n]=size(A);
if m ~= n || isequal(A,A')==false 
    error('Input error! Input must be a symmetric matrix!')
end

P=eye(n);
L_temp=zeros(n);
L=eye(n);
T=A;
for i=1:n-2
    r=maxIndex(T(:,i),i+1,n); 
    if abs(T(r,i))>eps
        if r~=i+1
            T([i+1,r],:)=T([r,i+1],:); % row (i+1) <-> row r
            T(:,[i+1,r])=T(:,[r,i+1]); % column (i+1) <-> column r
            I=eye(n);I(r,r)=0;I(i+1,r)=1;I(r,i+1)=1;I(i+1,i+1)=0;% construct P(i+1,r)
            P=I*P; %   P=P(n-2)...P(2)P(1)
        end
    
    end

%      _ _
% A <= GAG' ... until A is tridiagonal 


    e=zeros(1,n); e(i+1)=1; % unit row vector e
    L_temp(i+2:n,i)=T(i+2:n,i)/T(i+1,i);
    G=eye(n)-L_temp(:,i)*e;
    g=eye(n)+L_temp(:,i)*e; % g=inv(G)
    T=G*T*G.';
    L=L*g;
% It is worthing noting that L(1)(g(1)) is not equal to eye(n)+ X,
% where X(:,1)=[0;0;T(3:n,1)]/T(2,1)] and rest of X is 0, but equal
% to eye(n)+X, where X(:,2)=[0;0;T(3:n,1)]/T(2,1)] and rest is 0.
% Indeed, debug this source file, we can identify this subtle trap.
% Overall, the columns of L are shifted left a single unit in Aasen0.m!!! 
% And that's the very bug that "Aasen0.m" exists



end


end
