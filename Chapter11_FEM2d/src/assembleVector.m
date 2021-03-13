function phi = assembleVector(p,t,f)
nCells=size(t,1);
phi=zeros(size(p,1),1);
for k=1:nCells
    E=t(k,:); % element by index
    phi(E)=phi(E)+getElementVector(p(E,:),f);
end
end

function phi_loc = getElementVector(triangle,f)
phi_loc=f(triangle(:,1),triangle(:,2));
phi_loc=getArea(triangle)/3*phi_loc;
end

