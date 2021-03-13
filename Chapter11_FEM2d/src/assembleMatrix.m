% $ -\Delta u + ku = f $, @para mass_coef = $k$
% $ \int_{\Omega}\, kuv\;dx $ is termed as mass matrix
function A = assembleMatrix(p, t, mass_coef)
if nargin==2
    mass_coef=0;
end
N=size(p,1);
nCells=size(t,1);
ii=zeros(9*nCells,1); % prealloacate memory
jj=ii;
vv=ii;
idx=0;
for k=1:nCells
    E=t(k,:); % triangle element by index
    Ak=getElementMatrix(p(E,:), mass_coef);
    for i=1:3
        for j=1:3
            idx=idx+1;
            ii(idx)=E(i); jj(idx)=E(j);
            vv(idx)=Ak(i,j);
        end
    end
end
A=sparse(ii,jj,vv,N,N);
end

function Ak = getElementMatrix(triangle, mass_coef)
Ak=ones(3);
Ak(:,2:3)=triangle;
Ak=Ak\eye(3); % inv(Ak)
Ak=Ak(2:3,:);
area=getArea(triangle);
Ak=area*(Ak.')*Ak;
if mass_coef ~= 0
    Ak=Ak+mass_coef*area/12*[2 1 1;1 2 1; 1 1 2];
end
end

