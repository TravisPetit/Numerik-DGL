clc; clear all;

phi = @(x,y) -(x*y + x^2 + y^2);
c = -1;
xs = -linspace(-3,3,250);
ys = -linspace(-3,3,250);
p = 50;
[C, H, G, B] = makeGrid(phi, c, xs, ys, p);

figure();
scatter(G(1,:), G(2,:), "x");
hold on
scatter(B(1,:), B(2,:));
hold off
