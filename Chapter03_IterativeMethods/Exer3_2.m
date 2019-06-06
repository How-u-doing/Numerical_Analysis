% Exer 3.2
% 10170437 Mark Taylor
function Exer3_2()
    A=[1, 0.9, 0.9; 0.9, 1, 0.9; 0.9, 0.9, 1];  
    x=[4,3,7].';
    b=A*x; 
    w=[-1,-0.5,0.1,0.3,0.5,  0.6897,  0.71,0.72,0.8,1.2];
    %                           ^
    %                           |
    %                        w_opt=2/(0.1+2.8)=0.6897
    tol=1.0e-9;
    N=2000;
    X0=zeros(3,1); 
    for i=1:length(w)
        fprintf('Test(%d): w=%.4f\n',i,w(i));
        [x, k]=Richardson(A, b, w(i), tol, N, X0)
    end
    
    fprintf('Pick up a small neighborhood centered at 0.6897(posssible optimum omega)\n');
    fprintf('with an interval of 0.005 and the respective results are as followed:\n');
    w=0.6797:0.0050:0.6997;
    for i=1:length(w)
        fprintf('\nw=%.4f\n',w(i));
        [x, k]=Richardson(A, b, w(i), tol, N, X0)
    end 
    
    
    %###########################################################################
    % I got the w_opt=0.6825 which is not equal to its theoretical value 0.6897
    % may due to the accumulation of rounding error. 
    %###########################################################################
    fprintf('Now let''s choose another neighborhood centered at 0.6825 with \n');
    fprintf('an interval of 0.0005 and the respective results are as followed:\n');
    w=0.6810:0.0005:0.6840;
    for i=1:length(w)
        fprintf('\nw=%.4f\n',w(i));
        [x, k]=Richardson(A, b, w(i), tol, N, X0)
    end
       
end