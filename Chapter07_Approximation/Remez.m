% Remez's 1st algorithm£¬to find an optimal polynomial that can best approx f(x)
% 10170437 Mark Taylor
function [pn,Enf,x,k]=Remez(f,x,y,a,b,epsilon,maxIt)
%INPUT:
    % f: the given function supporting vector input, e.g. f=@(x)sin(x)
    % x: a vecter containing n+2 initial test points
    % y: corresponding function values of x, i.e. f(x)
    % [a,b]: the given interval
    % epsilon: the given tolerance, stopping criterion
    % maxIt: maximum iteration
%OUTPUT:
    % pn: optimal polynomial
    % Enf: maximum error
    % x: determined alternative points
    % k: no. of iterations
    
% set default input arguments
if nargin<7
    maxIt=100;
    if nargin<6
        epsilon=1e-6;
        if nargin<5
            error('Error! Insufficient input arguments!\n')
        end
    end
end
if isrow(x)
    x=x.';
end
if isrow(y)
    y=y.';
end
n=length(x)-2; % the optimal polynomial to f(x) of degree at most n
A=zeros(n+2);  % coefficient matrix
A(:,2)=ones(n+2,1);     % coefficient of a0 (=1)
A(1,1)=1;
for i=2:(n+2)    
    A(i,1)=-1*A(i-1,1); % alternative sign for sigma (En(f) or -En(f))  
    % or A(i,1)=1-2*i%2 (mod(i,2)) , i=1,2,...,n+2
end
for j=1:n
    A(:,j+2)=x.^j;    
end
% X0=[sigma, a0, a1, ... , an].' 
X0=A\y;        % initial solutions

syms t;
s=1;pn=0;
% construct initial pn(t)
for i=0:n    
    pn=pn+X0(i+2)*s;
    s=s*t;
end

for k=1:maxIt
    g=pn-f(t);  % error function g(t)=pn(t)-f(t)
    
    % compute sup g(t), a<=t<=b, and find xstar which is the sup point
    xslv=solve(diff(g,t)==0,t); % all solutions (may lie in outside [a,b])
    xslv=double(xslv);
    n0=length(xslv);
    ga=abs(subs(g,t,a));
    gb=abs(subs(g,t,b));
    if ga>gb
        sup=ga;
        xstar=a;
    else
        sup=gb;
        xstar=b;
    end
    for i=1:n0
        if isreal(xslv(i)) && a<xslv(i) && xslv(i)<b
            gi=abs(subs(g,t,xslv(i)));
            if gi>sup
                sup=gi;
                xstar=xslv(i);
            end
        else % not in (a,b)
            continue;
        end    
    end
    
    % substitute with xstar, update x & A
    gxs=subs(g,t,xstar);
    for i=1:n+2
        if abs(xstar-x(i))<=epsilon
            Enf=vpa(sup,7);
            pn=vpa(pn,7);
            return;
        end
    end
    for i=1:n+1
        % case 1: most likely to happen
        if x(i)<xstar && xstar< x(i+1)
            if gxs*subs(g,t,x(i))>0
                x(i)=xstar;
                y(i)=f(xstar);
                for j=1:n
                    A(i,j+2)=xstar^j;
                end                
            else
                x(i+1)=xstar;
                y(i+1)=f(xstar);
                for j=1:n
                    A(i+1,j+2)=xstar^j;
                end
            end
        end 
    end
    % case 2
    if a<= xstar && xstar<x(1)
        if gxs*subs(g,t,x(1))>0
            x(1)=xstar;
            y(1)=f(xstar);
            for j=1:n
                A(1,j+2)=xstar^j;
            end
        else            
            x(2:n+2)=x(1:n+1);
            x(1)=xstar;
            y=f(x);
            for j=1:n
                A(:,j+2)=x.^j;    
            end
        end                
    end
    % case 3
    if x(n+2)<xstar && xstar<=b
        if gxs*subs(g,t,x(n+2))>0
            x(n+2)=xstar;
            y(n+2)=f(xstar);
            for j=1:n
                A(n+2,j+2)=xstar^j;
            end
        else            
            x(1:n+1)=x(2:n+2);
            x(n+2)=xstar;
            y=f(x);
            for j=1:n
                A(:,j+2)=x.^j;    
            end
        end                
    end
        
    X1=A\y;    
    % update pn(t)
    s=1;pn=0;
    for i=0:n    
        pn=pn+X1(i+2)*s;
        s=s*t;
    end
    
    % tell when to stop
    stop=true;
    for i=0:n
        if abs(X1(i+2)- X0(i+2))>epsilon
            stop=false;
            break;
        end
    end
    if stop==true
        Enf=vpa(sup,7);
        pn=vpa(pn,7);
        return;
    else % update X0
        X0=X1;
    end
end
fprintf('Maximum iteration exceeded£¡The results in the last iteration as followed:\n')
Enf=vpa(sup,7);
pn=vpa(pn,7);
end

