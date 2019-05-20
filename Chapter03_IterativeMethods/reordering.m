% Interchange rows to let digonal elements nonzero(large enough)
% 10170437 Mark Taylor
function [isSuccessful, isReordered, newA]=reordering(A)
m=size(A,1);
isSuccessful=true;
isReordered=false;
for i=1:m
% Interchange rows when diagonal elements are too small
    if abs(A(i,i))<=1.0e-6
        isSuccessful=false;
        for j=1:m
             if abs(A(j,i))>1.0e-6 && abs(A(i,j))>1.0e-6
                 A([i,j],:)=A([j,i],:);
                 isSuccessful=true;
                 isReordered=true;
                 break;
             end
        end
    end        
end   

newA=A;
end


