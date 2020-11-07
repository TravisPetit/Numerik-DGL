%% explizitRK
% Solve the ODE y'(x) = f(x, y(x)) using the explicit Runge-Kutta method where
% y: [a,b] -> R^d
%
% Parameters:
% -FUNCTION-         f: [a,b] x R^d -> R^d
% -FLOAT-            a  y and the first component of f are defined over [a,b]
% -FLOAT-            b @see a
% -dx1 FLOAT-        ya  the initial value y(a)
% -POS INT-          n  number of estimated values for y in the interval [a,b]
% -mx1 FLOAT-        alpha  the values alpha1, ..., alpham where m is the consitency order
% -m-1 x m-1 FLOAT-  beta  the values beta_ij for i = 2,...,m and j = 1,...,m-1 (beta_1: and beta_m: are not needed because they are 0)
% -mx1 FLOAT-        gamma  the values gamma_1, ..., gamma_m
%
% Output:
% -1 x n+1 FLOAT- x  evenly spaced points between a and b
% -d x n+1 FLOAT- y  = [eta_0, ..., eta_n] the estimations of [y(x_0), ... , y(x_n)]
function [x, y] = explizitRK(f, a, b, ya, n, alpha, beta, gamma, h)
    x = linspace(a, b, n+1);
    y = ya;
    [m, ~] = size(beta);
     
    for i = 1:n
        % "k" = k_1^i
        k = f( [x(i) + h* alpha(1) ; y(:,i)] );
        for j = 2:m
            % linear combination beta_j,1 k_1 + ... + beta_j,j-1 * k_j-1
            % written as
            % [ k1 ... kj-1 ] (in R^m x j-1)  * [beta_j,1 ..., beta_j,j-1] (in R^j-1 x 1)
            sum_beta_k = k(:, 1:j-1) * beta(j, 1:j-1)';
          
            inner = [x(i) + h * alpha(j); y(:,i) + h * sum_beta_k];
            k(:,j) = f(inner);
        end
         y(:, i+1) = y(:,i) + h* k * gamma';
    end
end      
                
                
                
                