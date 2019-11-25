% Composite Integration 
function T=CompositeInt(f,a,b,n,method)

h=(b-a)/n;
T=0; S=f(a+1/2*h);
for i=1:n-1
    T=T+f(a+i*h);
    S=S+f(a+(i+1/2)*h);
end
T=f(a)+2*T+f(b);
if method=="trapezium"
    T=h/2*T;
    return;
elseif method=="Simpson"
    T=h/6*(4*S+T);
    return;     
else
    erorr('Method must be either "trapezium" or "Simpson"!')
end
end

