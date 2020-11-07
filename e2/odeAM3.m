%% odeAB2
% Solve the ODE y'(x) = f(x, y(x)) using the 3-step Adams-Moulton method where
% f: [a,b] x R^d -> R^d  and  y: [a,b] -> R^d
%
% Parameters:
% -FUNCTION-   f: [a,b] x R^d -> R^d
% -FUNCTION-   fy: [a,b] x R^d -> R^(d x d)
% -FLOAT- a    y and the first component of f are defined over [a,b]
% -FLOAT- b    @see a
% -dx1 FLOAT-  ya  the initial value y(a)
% -POS INT-    n  number of estimated values for y in the interval [a,b
% -FLOAT-      sigma  used in the newton method that is used to solve non linear equations on each iteration
% -FLOAT-      tol  see @sigma
% -POS INT-    maxiter  see @sigma
%
% Output:
% -1 x n+1 FLOAT- x  evenly spaced points between a and b
% -d x n+1 FLOAT- y  = [eta_0, ..., eta_n] the estimations of [y(x_0), ... , y(x_n)]
function [x, y] = odeAM3(f, fy, a, b, ya, n, sigma, tol, maxiter)
    x = linspace(a, b, n+1);
    h = (b-a)/n;
    theta = 1/2;

    % bootstrap
    [~, y] = thetaSchema(f, fy, a, b, ya, 2, theta, sigma, tol, maxiter, h);
    
    [d, ~] = size(y);
    I = eye(d);
    for i = 1:n-2
        t1 = 1/24 * f( [x(i); y(:,i)] );
        t2 = -5/24 * f( [x(i+1); y(:,i+1)] );
        t3 = 19/24 * f( [x(i+2); y(:,i+2)] );
        
        t4 = @(z) 3/8 * f( [x(i+3); z] );
        t4_prime = @(z) 3/8 * fy( [x(i+3); z] );
        
        phi = @(z) y(:, i+2) + h * (t1 + t2 + t3 + t4(z)) - z;
        phi_prime = @(z) h * t4_prime(z) - I;
        
        y(:, i+3) = newtonIterationSWS(phi, phi_prime, y(:,i+1), sigma, tol, maxiter);
    end
end