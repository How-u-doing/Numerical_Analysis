function phi = assembleVectorByGaussQuad(p, t, f)
nCells = size(t,1);
phi = zeros(size(p,1),1);
for k = 1:nCells
    vidx = t(k,:); % global vertex indices of k-th cell
    phi(vidx) = phi(vidx)+getElementVector(p,vidx,f);
end
end

function phi_loc = getElementVector(p, vidx, f)
% 6-point quadrature rule of order 4 on the unit triangle [0 0;1 0;0 1]
w = [1/60; 1/60; 1/60; 9/60; 9/60 ; 9/60];
x = [1/2 0; 1/2 1/2; 0 1/2; 1/6 1/6; 1/6 2/3; 2/3 1/6];
ref_shape_val = [1-x(:,1)-x(:,2), x(:,1), x(:,2)]; % values on quadrature points
    
bK = p(vidx(1),:);
BK = [p(vidx(2),:)-bK; p(vidx(3),:)-bK];
det_BK = abs(det(BK));

% transform quadrature points
y = x*BK+bK;

f_val = f(y(:,1),y(:,2));
phi_loc = sum(w.*f_val.*ref_shape_val)*det_BK;
phi_loc = phi_loc.';
end

