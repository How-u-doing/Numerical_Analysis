% Hermite Iterpolation test

fprintf('test for x^2 :')
X=[-1,0,1],Y=[1,0,1],Yd=[-2,0,2]
[y0,f] = HI(X,Y,Yd,[3,4,5])

fprintf('\n\nconstruct a polynomial that is tangent to sin(x) at 0 & pi (the answer is H(x) = x - x^2/pi) :')
X2=[0,pi],Y2=[0,0],Yd2=[1,-1]
[y1,h] = HI(X2,Y2,Yd2,pi/2)
[H, newTable] = Hermite([X2.',Y2.',Yd2.'])
