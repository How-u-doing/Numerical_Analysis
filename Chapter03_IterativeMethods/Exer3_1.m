% Exer 3.1
% 10170437 Mark Taylor
function Exer3_1()
A1=[1, 0.9, 0.9; 0.9, 1, 0.9; 0.9, 0.9, 1];
A2=[1, 0, 1; -1, 1, 0; 1, 2, -3];
x=[4,3,7].';
b1=A1*x;
b2=A2*x;
s=[0,3,10,100];
tol=1.0e-9;
N=2000;
for i=1:4
    fprintf('*******************************************\n');
    fprintf('*                Part(%d)                  *\n',i)
    fprintf('*******************************************\n');
    x_0=[s(i),s(i),s(i)].'  
    fprintf('Test for A1:\n')  
    fprintf('Results of Jacobi Iteration:\n');
    [y1J, k1J]=Jacobi_C(A1, b1, tol, N, x_0)
    fprintf('Results of Gauss-Seidel Iteration:\n');
    [y1G, k1G]=GaussSeidel_C(A1, b1, tol, N, x_0)
    
    fprintf('################################\n');
    fprintf('Test for A2:\n')  
    fprintf('Results of Jacobi Iteration:\n');
    [y2J, k2J]=Jacobi_C(A2, b2, tol, N, x_0)
    fprintf('Results of Gauss-Seidel Iteration:\n');
    [y2G, k2G]=GaussSeidel_C(A2, b2, tol, N, x_0)
end

end
