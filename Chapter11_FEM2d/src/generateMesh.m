function [p,t]=generateMesh(x,y,n_x,n_y)
h_x=(x(2)-x(1))/n_x;
h_y=(y(2)-y(1))/n_y;
[X,Y] = meshgrid(x(1):h_x:x(2), y(1):h_y:y(2));
N=(n_x+1)*(n_y+1);
x=reshape(X,N,1);
y=reshape(Y,N,1);
p=[x,y];
t=delaunay(x,y);
end