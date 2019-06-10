% Givens Matrix 
% 10170437 Mark Taylor
function G = givens_cs(x_i, x_j)
%     _              _
%  G=|   cos    sin   |
%    |_ -sin    cos  _|

sq=sqrt(x_i^2+x_j^2);
c=x_i/sq;
s=x_j/sq;
G=eye(2);
G(1,1)=c;
G(1,2)=s;
G(2,1)=-s;
G(2,2)=c;
end
