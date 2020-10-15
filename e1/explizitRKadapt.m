%% explizitRKadapt
% Solve the ODE y'(x) = f(x, y(x)) using the explicit Runge-Kutta method where
% y: [a,b] -> R^d
% with an adaptative step size.
%
% Parameters:
% -FUNCTION-         f: [a,b] x R^d -> R^d
% -FLOAT-            a  y and the first component of f are defined over [a,b]
% -FLOAT-            b @see a
% -dx1 FLOAT-        ya  the initial value y(a)
% -FLOAT-            h0  the initial step size
% -POS INT-          n  number of estimated values for y in the interval [a,b]
% -mx1 FLOAT-        alpha  the values alpha1, ..., alpham where m is the consitency order
% -m-1 x m-1 FLOAT-  beta  the values beta_ij for i = 2,...,m and j = 1,...,m-1 (beta_1: and beta_m: are not needed because they are 0)
% -mx1 FLOAT-        gamma  the values gamma_1, ..., gamma_m
% -mx1 FLOAT-        gammah  the values gammahat_1, ..., gammahat_1 used in the control procedure
% -POS INT-          p  consistency order of the vanilla RK procedure
% -FLOAT-            tau  constant factor < 1 by which the control procedure outperforms the tolerance-over-approximation ratio (I think)
function [x, y] = explizitRKadapt(f, a, b, ya, h0, alpha, beta, gamma, gammah, p, tau, epsilon)
    y = ya;
    x = a;
    h_i = h0;
    delta_i = epsilon;
    
    [m, ~] = size(beta);
    exponent = 1/(p+1);
    
    i = 1;
    
    while(x(i) < b)
        while true
            h_i = tau * (epsilon/delta_i)^exponent * h_i;
            
            % don't allow for step sizes that would go beyond the bounds
            h_i = min(h_i, b - x(i));
            
            k_i = evaluate(f, x(i), h_i, alpha, y(:,i), beta, m);
            y(:,i+1) = y(:,i) + h_i * k_i * gamma';
            delta_i = h_i * norm(k_i * (gamma - gammah)', 2);
            
            if delta_i <= epsilon
                break
            end
        end
        x(i+1) = x(i) + h_i;
        i = i + 1;
    end
    
end

% compute the K matrix.
function k_i = evaluate(f, x_i, h_i, alpha, eta_i, beta, m)
    k_i = f( [x_i + h_i * alpha(1) ; eta_i] );
    for j = 2:m
        sum_beta_k = k_i(:, 1:j-1) * beta(j, 1:j-1)';      
        inner = [x_i + h_i * alpha(j); eta_i + h_i * sum_beta_k];
        k_i(:,j) = f(inner);
    end
end