% Givens Matrix 
% 10170437 Mark Taylor
function T = Givens(x_i, x_j,i,j,n)
sq=sqrt(x_i^2+x_j^2);
%                  Givens Rotation
%
c=x_i/sq;%     / 1                  \
s=x_j/sq;%     |   1                |
T=eye(n);%     |     c     s        |   <-- row i
%              |                    |
%        T  =  |                    |
%              |    -s     c        |   <-- row j
%              |              1     |
T(i,i)=c;%     |                 1  |
T(i,j)=s;%     \                    /
T(j,i)=-s;%         ¡ü     ¡ü
T(j,j)=c;%    column i     j
end
