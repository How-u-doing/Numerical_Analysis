% Givens Rotation Matrix 
% 10170437 Mark Taylor
function P = givens(x_i, x_j)
%     _              _      
%  P=|   cos    sin   |    
%    |_ -sin    cos  _|

sq=sqrt(x_i^2+x_j^2);
c=x_i/sq;
s=x_j/sq;
P=eye(2);
P(1,1)=c;
P(1,2)=s;
P(2,1)=-s;
P(2,2)=c;
end
