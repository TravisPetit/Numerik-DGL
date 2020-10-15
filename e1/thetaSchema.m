%% thetaSchema
% Solve the ODE y'(x) = f(x, y(x)) with the theta schema where
% y: [a,b] -> R^d
%
% Parameters
% -FUNCTION-   f: [a,b] x R^d -> R^d
% -FUNCTION-   fy: [a,b] x R^d -> R^(d x d)  derivative of f in the second variable
% -FLOAT- a    y and the first component of f are defined over [a,b]
% -FLOAT- b    see @a
% -dx1 FLOAT-  ya  the initial value y(a)
% -POS INT-    n  number of estimated values for y in the interval [a,b]
% -FLOAT-      theta  
% -FLOAT-      sigma  used in the newton method that is used to solve non linear equations on each iteration
% -FLOAT-      tol  see @sigma
% -POS INT-    maxiter  see @sigma
%
% Output:
% -1 x n+1 FLOAT- x = [x_0, ... ,x_n] evenly distributed between a and b.
% -d x n+1 FLOAT- y = [eta_0, ..., eta_n] the estimations of [y(x_0), ... , y(x_n)]
function [x, y] = thetaSchema(f, fy, a, b, ya, n, theta, sigma, tol, maxiter)
    x = linspace(a, b, n+1);
    h = (b-a)/n;
    y = ya;
    [d, ~] = size(y);
    I = eye(d);
    
    for i = 1:n
        % phi: R^d -> R^d
        phi = @(q) q - y(:,i) - h * ( (1-theta) * f( [x(i); y(:,i)] ) + theta * f( [x(i+1); q] ) );
        
        %phi_prime: R^d -> R^d x R^d
        phi_prime = @(q) I - h * theta * fy( [x(i+1); q] );
        
        y(:, i+1) = newtonIterationSWS(phi, phi_prime, y(:,i), sigma, tol, maxiter);
    end
end





