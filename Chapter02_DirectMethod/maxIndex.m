function index=maxIndex(x,first,last)
% Return the index of component with the largest absolute value that 
% ranges from first to last in vector x
[m,n]=size(x);
if min(m,n)~=1
    error('First input argument must be a vector!');
end

if first>=1 && last<=length(x) && first<=last
    MAX=abs(x(first));
    index=first;
    for k=first+1:last
        if abs(x(k))>MAX
            MAX=abs(x(k));
            index=k;
        end
    end
else
    error('first and last must satisfy:first>=1&&last<=length(x)&&first<=last');
end
end

