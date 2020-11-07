%% odeAB2
% Solve the ODE y'(x) = f(x, y(x)) using the 3-step Adams-Bashforth method where
% f: [a,b] x R^d -> R^d  and  y: [a,b] -> R^d
%
% Parameters:
% -FUNCTION-   f: [a,b] x R^d -> R^d
% -FLOAT- a    y and the first component of f are defined over [a,b]
% -FLOAT- b    @see a
% -dx1 FLOAT-  ya  the initial value y(a)
% -POS INT-    n  number of estimated values for y in the interval [a,b]
%
% Output:
% -1 x n+1 FLOAT- x  evenly spaced points between a and b
% -d x n+1 FLOAT- y  = [eta_0, ..., eta_n] the estimations of [y(x_0), ... , y(x_n)]
function [x, y] = odeAB3(f, a, b, ya, n)
    x = linspace(a, b, n+1);
    h = (b-a)/n;
    
    RK_alpha = [0, 1/2, 1/2, 1];
    RK_beta = [0,   0,   0, 0;
            1/2, 0,   0, 0;
            0,   1/2, 0, 0;
            0,   0,   1, 0];
    RK_gamma = [1/6, 1/3, 1/3, 1/6];
    
    [~, y] = explizitRK(f, a, b, ya, 2, RK_alpha, RK_beta, RK_gamma, h);
    
    tmp = [f( [x(1);y(:,1)] ), f( [x(2); y(:,2)] ), f( [x(3); y(:,3)] )];
    
    for i = 1:n-2
        y(:,i+3) = y(:,i+2) + h * (5/12 * tmp(:,1) - 4/3 * tmp(:,2) + 23/12 * tmp(:,3));
        tmp = [tmp(:,2), tmp(:,3), f( [x(i+3); y(:,i+3)] )];
    end
end