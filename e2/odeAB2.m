%% odeAB2
% Solve the ODE y'(x) = f(x, y(x)) using the 2-step Adams-Bashforth method where
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
function [x, y] = odeAB2(f, a, b, ya, n)
    x = linspace(a, b, n+1);
    h = (b-a)/n;
    
    RK_alpha = [0, 1/2, 1/2, 1];
    RK_beta = [0,   0,   0, 0;
            1/2, 0,   0, 0;
            0,   1/2, 0, 0;
            0,   0,   1, 0];
    RK_gamma = [1/6, 1/3, 1/3, 1/6];
    
    % bootstrap
    [~, y] = explizitRK(f, a, b, ya, 1, RK_alpha, RK_beta, RK_gamma, h);
    
    tmp = [f( [x(1);y(:,1)] ), f( [x(2); y(:,2)] )];
    
    for i = 1:n-1
        %y(:,i+2) = y(:,i+1) + h * (-1/2 * f( [x(i); y(:,i)] ) + 3/2 * f( [x(i+1); y(:,i+1)] ));
        y(:,i+2) = y(:,i+1) + h * (-1/2 * tmp(:,1) + 3/2 * tmp(:,2));
        tmp = [tmp(:,2), f( [x(i+2); y(:,i+2)] )];
    end
end