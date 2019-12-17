% Composite Integration 
function T=CompositeInt(f,a,b,n,method,m)

h=(b-a)/n;
T=0;
for i=1:n-1
    T=T+f(a+i*h);
end
T=f(a)+2*T+f(b);
if method=="trapezoid"
    T=h/2*T;
    return;
elseif method=="Simpson"
    S=0;
    for i=0:n-1
        S=S+f(a+(i+1/2)*h);         
    end
    T=h/6*(4*S+T);
    return;     
elseif method=="Gauss"
    T=0;
    x=a:h:b;
    for i=1:n
        T=T+Gauss_quadrature(f,x(i),x(i+1),m);
    end
    return
else
    erorr('Method must be either "trapezium" or "Simpson" or "Gauss"!')
end
end

