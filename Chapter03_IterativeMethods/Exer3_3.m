% Exer 3.3
% 10170437 Mark Taylor
function Exer3_3()
C = 5*eye(6) + rand(6);
A = C.' * C;
x=[1 7 0 4 3 7].';
b=A*x;
tol=1.0e-9;
N=1000;
fprintf('Gradient Descent method\n');
[yGD, rGD, kGD]=G_D(A, b, tol, N)

fprintf('Minimum Residual method\n');
[yMR, rMR, kMR]=M_R(A, b, tol, N)

fprintf('Conjugate Gradient method\n');
[yCG, rCG, kCG]=C_G(A, b, tol, N)

end
