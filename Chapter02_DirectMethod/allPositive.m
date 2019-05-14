function tf= allPositive(x)
% Tell if each component of vector x is positive 
tf=true;
for i=1:length(x)
    if x(i)<=eps % Generally speaking, eps=2.2204e-016 >0
        tf=false;
        break;
    end
    
end

end
