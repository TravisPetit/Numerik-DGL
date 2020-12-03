%% makeGrid
% Creates a raster for the level set defined by phi and c.
%
% Parameters:
% -FUNCTION- phi: R^2 -> R, used to define the level set {(x,y) in R^2 : phi(x,y) < c}
% -FLOAT- c see @phi
% -1 x r FLOAT- xs partition of the x axis
% -1 x s FLOAT- ys partition of the y axis
% -POS INT- p number of bisection steps
%
% Output:
% -2 x nd FLOAT- G column vectors = inner raster points
% -2 x nb FLOAT- B column vectors = the boundary points
% -4 x nd INT- C column vectors = index of left, right, upper and lower point
% positive => inner point => index of G
% negative => boundary point => index of B
% -4 x nd FLOAT- H matrix that has the same structure as C but instead of the indices it contains the distances
function [C, H, G, B] = makeGrid(phi, c, xs, ys, p)
    [~, r] = size(xs);
    [~, s] = size(ys);
    I = zeros(r,s);

    % compute I
    nd = 0;
    for j = 1:r
        for k = 1:s
            if phi(xs(j), ys(k)) > c
                nd = nd + 1;
                I(j,k) = nd;
            end
        end
    end
    
    % compute C part 1/2
    % and compute G
    C = zeros(4, nd);
    G = zeros(2, nd);
    for j = 2:r-1
        for k = 2:s-1
            if I(j,k) ~= 0
                C(1,I(j,k)) = I(j-1,k);
                C(2,I(j,k)) = I(j+1,k);
                C(3,I(j,k)) = I(j-1,k-1);
                C(4,I(j,k)) = I(j-1,k+1);
                
                G(1, I(j,k)) = xs(j);
                G(2, I(j,k)) = ys(k);
            end
        end
    end
        
    % compute C part 2/2 and B, H
    nb = nnz(~C); % number of zeros in C
    B = zeros(2, nb);
    m = 0;
    for j = 2:r-1
        for k = 2:s-1
            if I(j,k) ~=0
                
                if I(j-1,k) == 0
                    m = m + 1;
                    C(1, I(j, k)) = -m;
                    
                    left = xs(j);
                    right = xs(j-1);
                    for q = 1:p
                        mid = (right - left)/2;
                        if phi(mid, ys(k)) > c
                            right = mid;
                        else
                            left = mid;
                        end
                    end
                    B(1,m) = left;
                    B(2,m) = ys(k);
                end
                
                if I(j+1,k) == 0
                    m = m + 1;
                    C(2, I(j, k)) = -m;
                    
                    left = xs(j+1);
                    right = xs(j);
                    for q = 1:p
                        mid = (right - left)/2;
                        if phi(mid, ys(k)) <= c
                            right = mid;
                        else
                            left = mid;
                        end
                    end
                    B(1,m) = right;
                    B(2,m) = ys(k);
                end
                
                if I(j,k-1) == 0
                    m = m + 1;
                    C(3, I(j, k)) = -m;
                    
                    up = ys(k);
                    down = ys(k-1);
                    for q = 1:p
                        mid = (up - down)/2;
                        if phi(xs(j), mid) > c
                            down = mid;
                        else
                            up = mid;
                        end
                    end
                    B(1,m) = xs(j);
                    B(2,m) = up;
                end
                
                if I(j,k+1) == 0
                    m = m + 1;
                    C(4, I(j, k)) = -m;
                    
                    up = ys(k+1);
                    down = ys(k);
                    for q = 1:p
                        mid = (up - down)/2;
                        if phi(xs(j), mid) > c
                            up = mid;
                        else
                            down = mid;
                        end
                    end
                    B(1,m) = xs(j);
                    B(2,m) = down;
                end
                
            end
        end
    end
    
    H = zeros(4, nd);
end