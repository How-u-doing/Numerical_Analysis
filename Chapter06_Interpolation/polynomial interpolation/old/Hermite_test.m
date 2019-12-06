% Hermite test
% 10170437 Mark Taylor

table=[1.3,0.6200860,-0.5220232;1.6,0.4554022,-0.5698959;1.9,0.2818186,-0.5811571]
[H, newTable] = Hermite(table)
fprintf('H£¨1.5£©=\n');
vpa(subs(H,1.5),7)
