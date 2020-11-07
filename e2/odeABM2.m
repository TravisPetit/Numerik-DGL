%% odeAM2
% Solve the ODE y'(x) = f(x, y(x)) using the 2-step Adams-Bashforth method with a 2-step Adams-Moulton correction method where
% f: [a,b] x R^d -> R^d  and  y: [a,b] -> R^d
%
% Parameters:
% -FUNCTION-   f: [a,b] x R^d -> R^d
% -FLOAT- a    y and the first component of f are defined over [a,b]
% -FLOAT- b    @see a
% -dx1 FLOAT-  ya  the initial value y(a)
% -POS INT-    n  number of estimated values for y in the interval [a,b]
% -POS INT-    m  number of corrector iterations
%
% Output:
% -1 x n+1 FLOAT- x  evenly spaced points between a and b
% -d x n+1 FLOAT- y  = [eta_0, ..., eta_n] the estimations of [y(x_0), ... , y(x_n)]
function [x, y] = odeABM2(f, a, b, ya, n, m)
    h = (b-a)/n;
    
    % predictor
    [x,y] = odeAB2(f, a, b, ya, n);
    
    tmp = [f( [x(1); y(:,1)] ), f( [x(2); y(:,2)] )];
    
    for i = 1:n-1
        t1 = -1/12 * tmp(:,1);
        t2 = 2/3 * tmp(:,2);
        
        for step = 1:m
            t3 = 5/12 * f( [x(i+2);y(:,i+2)] );
            
            % corrector
            y(:,i+2) = y(:,i+1) + h * (t1 + t2 + t3);
        end
        tmp = [tmp(:,2), f( [x(i+2); y(:,i+2)] )];
    end
end