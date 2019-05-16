% Inverse matrix of a nonsingular upper triangular matrix
% If input is a lower triangular matrix L, its inverse matrix inv_L=Inv_U(L.').' 

% 10170437 Mark Taylor
function V=Inv_U(U) % V*U=I
    if all(diag(U))==false
        error('U is a singular matrix which has no inverse matrix')
    end
    n=size(U,2);
    V=zeros(n);
    
    V(1,1)=1/U(1,1);
    for j=2:n
        V(j,j)=1/U(j,j);
        for i=1:j-1
            V(i,j)=-1/U(j,j)*V(i,1:j-1)*U(1:j-1,j);
        end
        
    end

end