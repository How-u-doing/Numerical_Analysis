function [L_inf_err,L2_err,H1_err] = errorEstimates(p,t,u,u_r,grad_u_r)

% 6-point quadrature rule of order 4 on unit triangle [0 0;0 1;1 0]
w = [1/60; 1/60; 1/60; 9/60; 9/60 ; 9/60];
x = [1/2 0; 1/2 1/2; 0 1/2; 1/6 1/6; 1/6 2/3; 2/3 1/6];

nCells = size(t,1);
unit_base = [1-x(:,1)-x(:,2), x(:,1), x(:,2)];
grad_unit_base = [-1 -1 1 0 0 1];

L_inf_err = max(abs(u_r(p(:,1),p(:,2))-u));
L2_err = 0; H1_err = 0;  
for k = 1 : nCells
    % We can also use `GaussTriaQuad()`, but the evaluation of the values
    % under those of points after the affine transformation is nontrivial,
    % since u_FE on K'th triangle is: u1*b1 + u2*b2 + u3*b3, where bi is a
    % plane based on the triangle and only the i'th point has height 1.
    % Hence, we adopt another approach which directly use the unit base
    % functions, the precedure is identical to `GaussTriaQuad()`.
    
    vidx = t(k,:);
    
    bK = p(vidx(1),:);
    BK = [p(vidx(2),:)-bK; p(vidx(3),:)-bK];
    det_BK = abs(det(BK));
    
    % transform quadrature points
    y = x*BK+bK;
    
    u_EX = u_r(y(:,1),y(:,2));
    grad_u_EX = grad_u_r(y(:,1),y(:,2));
    
    u_FE = u(vidx(1))*unit_base(:,1)+u(vidx(2))*unit_base(:,2)+ ...
                                     u(vidx(3))*unit_base(:,3);
    
    grad_u_FE = (u(vidx(1))*grad_unit_base(1:2)+ ...
                 u(vidx(2))*grad_unit_base(3:4)+ ...
                 u(vidx(3))*grad_unit_base(5:6))*(inv(BK)).';
    
    tmp = sum(w.*(u_EX-u_FE).^2)*det_BK;
    L2_err = L2_err+tmp;
    
    H1_err = H1_err+tmp+sum(w.*sum((grad_u_EX-grad_u_FE).^2,2))*det_BK;
    
end

L2_err = sqrt(L2_err);
H1_err = sqrt(H1_err);

end

