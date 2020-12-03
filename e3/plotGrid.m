function plotGrid(C, H, G, B)
    K = C;
    K(C < 0) = size(G, 2) - C(C < 0);
    P = [G, B];

    plot(G(1, :), G(2, :), 'ks', ...
         [G(1, :); P(1, K(1, :))], [G(2, :); P(2, K(1, :))], 'r-', ...
         [G(1, :); P(1, K(2, :))], [G(2, :); P(2, K(2, :))], 'k-.', ...
         [G(1, :); P(1, K(3, :))], [G(2, :); P(2, K(3, :))], 'r-', ...
         [G(1, :); P(1, K(4, :))], [G(2, :); P(2, K(4, :))], 'k-.', ...
         B(1, :), B(2, :), 'rx');
    axis('equal');
end