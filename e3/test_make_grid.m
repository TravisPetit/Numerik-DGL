clc; clear all;

%%

phi = @(x,y) -(x*y + x^2 + y^2);
c = -1;
xs = linspace(-1.5,1.5,20);
ys = linspace(-1.5,1.5,20);
p = 10;
[C, H, G, B] = makeGrid(phi, c, xs, ys, p);

figure();
scatter(G(1,:), G(2,:), "x");
hold on
scatter(B(1,:), B(2,:));
hold off


%%

% n = 2;
% [phi, c, xl, yl] = domains('ellipsering');
% xs = linspace(xl(1), xl(2), 2^(n+1)+1);
% ys = linspace(yl(1), yl(2), 2^(n+1)+1);
% p = 7+n;
% [C, H, G, B] = makeGrid(phi, c, xs, ys, p);
% 
% f = figure();
% plotGridMesh(C, H, G, B);
% 
% 
% 
% MAT = load("ringGrid_level_2.mat");
% B_M = MAT.B;
% C_M = MAT.C;
% G_M= MAT.G;
% H_M = MAT.H;
% 
% g = figure();
% plotGridMesh(C_M, H_M, G_M, B_M);

% B
% B_M