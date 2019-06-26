% Givens Rotation Matrix 
% 10170437 Mark Taylor
function T = GV(x_i, x_j,i,j,n)
sq=sqrt(x_i^2+x_j^2);
%            Givens Rotation Matrix form
%
c=x_i/sq; %     / 1                  \
s=x_j/sq; %     |   1                |
T=eye(n); %     |     c     s        |   <-- row i
%               |       1            |
%        T  =   |         1          |
%               |    -s     c        |   <-- row j
%               |             1      |
T(i,i)=c; %     |               1    |
T(i,j)=s; %     \                 1  /
T(j,i)=-s;%           ^     ^
T(j,j)=c; %           |     |
          %    column i     j
end
