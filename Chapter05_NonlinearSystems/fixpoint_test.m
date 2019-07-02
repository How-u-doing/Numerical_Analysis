% Test for "fixpoint.m"
% 10170437 Mark Taylor

function fixpoint_test()
% Here's an example:
%   Say we wanna solve x^2-3=0, <-> x=phi(x), set phi(x) following functions:
%   #1.    phi1(x)= x^2+x-3, 
%   #2.    phi2(x)= 3/x, 
%   #3.    phi3(x)= x-1/4*(x^2-3),
%   #4.    phi4(x)= 1/2*(x+3/x).
%
% Then we need to choose an appropriate initial value, say x0=2(Anyway, it 
% would be beneficial if we choose one close to exact solution(3^.5) enough.
% And you know there is a maximum delta, i.e. |x-3^.5|<delta, if x0 has a  
% distance to 3^.5 that is larger than the maximum delta, the whole iterative 
% procedure does NOT converge! 

x0=2;
phi1=@(x) x^2+x-3;
fprintf('/***phi1(x)\n')
[x1,k1]=fixpoint(phi1,x0)
fprintf('***/end of phi1(x)\n\n\n')

phi2=@(x) 3/x;
fprintf('/***phi2(x)\n')
[x2,k2]=fixpoint(phi2,x0)
fprintf('***/end of phi2(x)\n\n\n')

phi3=@(x) x-1/4*(x^2-3);
fprintf('/***phi3(x)\n')
[x3,k3]=fixpoint(phi3,x0)
fprintf('***/end of phi3(x)\n\n\n')

phi4=@(x) 1/2*(x+3/x);
fprintf('/***phi4(x)\n')
[x4,k4]=fixpoint(phi4,x0)
fprintf('***/end of phi4(x)\n\n\n')

% Now let's choose two different x0s such that one can get the answer, the 
% other's not.(|x1-3^.5|<delta, |x2-3^.5|>=delta)
fprintf('\nFor phi3(x)=x-1/4*(x^2-3), we set x1 & x2 the values 5.7 and 5.8 respectively.\n\n\n')
fprintf('For x1=5.7\n')
[x31,k31]=fixpoint(phi3,5.7)
fprintf('\n\nFor x2=5.8\n')
[x32,k32]=fixpoint(phi3,5.8)

end

