% $ -\nabla\cdot(\bm{\alpha}\nabla u) + \gamma u = f $
function A = assembleReactionDiffusionMatrix(p, t, alpha, gamma)
N=size(p,1);
nCells=size(t,1);
ii=zeros(9*nCells,1); % prealloacate memory
jj=ii;
vv=ii;
idx=0;
for k=1:nCells
    vidx=t(k,:); % global vertex indices of k-th cell
    Ak=getElementMatrix(p,vidx,alpha,gamma);
    for i=1:3        
        for j=1:3
            idx=idx+1;
            ii(idx)=vidx(i); jj(idx)=vidx(j);
            vv(idx)=Ak(i,j);
        end
    end
end
A=sparse(ii,jj,vv,N,N);
end

function Ak = getElementMatrix(p, vidx, alpha, gamma)
% 6-point quadrature rule of order 4 on the unit triangle [0 0;1 0;0 1]
w = [1/60; 1/60; 1/60; 9/60; 9/60 ; 9/60];
x = [1/2 0; 1/2 1/2; 0 1/2; 1/6 1/6; 1/6 2/3; 2/3 1/6];
ref_shape_val = [1-x(:,1)-x(:,2), x(:,1), x(:,2)]; % values on quadrature points
grad_ref_shape = [-1 1 0;  % gradients of reference shape functions
                  -1 0 1];
              
bK = p(vidx(1),:);
BK = [p(vidx(2),:)-bK; p(vidx(3),:)-bK];
det_BK = abs(det(BK));
inv_BK = BK\eye(2); % inv(BK)

% transform quadrature points
y = x*BK+bK;

gamma_val = gamma(y(:,1),y(:,2));
Ak = zeros(3);
for i=1:6 % loop over quadrature points
    alpha_val = alpha(y(i,1),y(i,2)); % a 2-by-2 symmetric matrix or a scalar
    
    trf_grad = inv_BK*grad_ref_shape; % transformed gradients, 2-by-3
    Ak_alpha = (alpha_val*trf_grad).'*trf_grad;
    
    Ak_gamma = gamma_val(i)*(ref_shape_val(i,:).')*ref_shape_val(i,:);
    
    Ak = Ak+w(i)*(Ak_alpha+Ak_gamma)*det_BK;
end
end
