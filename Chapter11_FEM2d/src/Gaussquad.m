% Gaussian quadrature
function T = Gaussquad(f,a,b,n)
if nargin<4
    n=8;
end
% Gaussian quarature points & weights on [-1,1]
if n == 2
    w = [1,1];                      % weights
    p = [-1/sqrt(3),1/sqrt(3)];     % points
elseif n == 4
    w = [0.3478548451,0.3478548451,0.6521451549,0.6521451549];
    p = [0.8611363116,-0.8611363116,0.3399810436,-0.3399810436];
elseif n == 8
    w = [0.1012285363,0.1012285363,0.2223810345,0.2223810345,0.3137066459,...
	         0.3137066459,0.3626837834,0.3626837834];
    p = [0.9602898565,-0.9602898565,0.7966664774,-0.7966664774,0.5255324099,...
			-0.5255324099,0.1834346425,-0.1834346425];
else
    error('n must be either 2 or 4 or 8')
end

% take a linear transform: x = t*(b-a)/2 + (a+b)/2, where t belongs to [-1,1]
% int(f,[a,b]) = (b-a)/2 * int(f(t*(b-a)/2 + (a+b)/2), [-1,1])
w = 0.5*(b-a)*w;
p = 0.5*(b-a)*p+0.5*(a+b);

T = w*f(p).';
end

