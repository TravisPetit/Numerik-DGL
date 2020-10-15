%% newtonIterationSWS
% Solve the nonlinear system of equations g(z) = 0
%
% Parameters:
% -FUNCTION-   g : R^d -> R^d
% -FUNCTION-   gp : R^d -> R^(d x d)  returns the jacobian for the given input
% -dx1 FLOAT-  z0  the initial guess
% -FLOAT-      sigma  in [0,1], the "strictness" of each step size. sigma = 0 is the classical approach.
% -FLOAT-      tol  the tolerance of error
% -POS INT-    maxiter   the maximum amount of iterations
%
% Output:
% z in R^d such that ||g(z)|| is close to 0.
function z = newtonIterationSWS(g, gp, z0, sigma, tol, maxiter)
    z = z0;
    iter = 0;
    while (iter < maxiter)
        d = gp(z) \ g(z);
        if (norm(d,2) <= tol) %i.e. if z_n+1 "=" z_n
            return
        end
        alpha = 1;
        while(norm(g(z - alpha*d), 2) > norm(g(z),2) * (1-sigma*alpha))
            alpha = alpha / 2;
        end
        z = z - alpha * d;
        iter = iter + 1;
    end
end
