clc; clear all;
%phi = @(x,y) -(x*y + x^2 + y^2);
% c = -1;
[phi, c, xl, yl] = domains('ellipsering');

xs = linspace(-1.5, 1.5, 500);
ys = linspace(-1.5, 1.5, 500);
p = 11;

g = @(x) 0;

[C, H, G, B] = makeGrid(phi, c, xs, ys, p);
[A, Ab] = shortleyWeller(C, H);

[~, num_points_boundary] = size(B);

[eigenvectors, lambda] = eigs(A, 10, "smallestabs");

v = zeros(num_points_boundary, 1);
[F, P] = makeMesh(C, H, G, B);

for i = 1:10
    subplot(2,5,i)
    plotMeshFunction(F, P, eigenvectors(:,i), v);
    title(sprintf("\\lambda_{%d} = %f", i, lambda(i,i)));
end