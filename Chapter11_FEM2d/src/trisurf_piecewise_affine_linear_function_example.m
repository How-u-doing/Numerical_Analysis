[x,y] = meshgrid((1:2:10)/10,(1:2:10)/10);
tri = delaunay(x,y);
z = peaks(5); 
z = z + .3*rand(5); % oscillation from the Gaussian distributions
zmax=max(max(z)); zmin=min(min(z));
a=.05; b =1;
z = a+(z-zmin)/(zmax-zmin)*(b-a); % transform to [a,b]

figure(1)
trisurf(tri,x,y,z)
figure(2)
triplot(tri,x,y)

% plot the mesh in 3d coordinate system
figure(1)
hold on
t=[.1, .9]; z0=[0,0];
for v=.1:.2:.9
    plot3(t,[v,v],z0,'r');
    plot3([v,v],t,z0,'r');
end
for v=.3:.2:.9
    plot3([.1,v],[v,.1],z0,'r')
end
for v=.3:.2:.7
    plot3([v,.9],[.9,v],z0,'r')
end

% You may want to adjust the figure to a satisfying position first
% fig2svg("../svg/trisurf_piecewise_affine_linear_function_example.svg")
%%% NOTE %%%
% When exporting to svg, the red mesh will cover (be on the upper layer) 
% part of the surface plot if we first draw the surface then the mesh.
% But we can use some svg editor (for example, Inkscape) to raise the layer
% of the surface so that we can use the svg in our paper.

