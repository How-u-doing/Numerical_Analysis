% Gaussian quadrature over a triangle
function S = GaussTriaQuad(triangle, f)
% @para triangle: a 3x2 matrix, 3 pairs of (x,y) coordinates
% @para f: function handle of a two-variale func to be integrated

% 6-point quadrature rule of order 4 on unit triangle [0 0;0 1;1 0]
w = [1/60; 1/60; 1/60; 9/60; 9/60 ; 9/60];
x = [1/2 0; 1/2 1/2; 0 1/2; 1/6 1/6; 1/6 2/3; 2/3 1/6];

% affine mapping to tranform unit triangle to a general triagle
%
%  |\        \Phi(x^) = F * x^ + \tau           /\
%  |  \     --------------------------->       /   \
%  |____\                                     /______\
%

bK = triangle(1,:);
BK = [triangle(2,:)-bK; triangle(3,:)-bK];
det_BK = abs(det(BK));

% transform quadrature points and integrate
x = x*BK+bK;
S = sum(w.*f(x(:,1),x(:,2)))*det_BK;
end

