%% odeAM2
% Solve the ODE y'(x) = f(x, y(x)) using the 3-step Adams-Bashforth method with a 3-step Adams-Moulton correction method where
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
function [x, y] = odeABM3(f, a, b, ya, n, m)
    h = (b-a)/n;
    
    % predictor
    [x,y] = odeAB3(f, a, b, ya, n);
    
    tmp = [f( [x(1);y(:,1)] ), f( [x(2); y(:,2)] ), f( [x(3); y(:,3)] )];
    
    for i = 1:n-2
        t1 = 1/24 * tmp(:,1);
        t2 = -5/24 * tmp(:,2);
        t3 = 19/24 * tmp(:,3);
        
        for step = 1:m
            t4 = 3/8 * f( [x(i+3);y(:,i+3)] );
            
            % corrector
            y(:,i+3) = y(:,i+2) + h * (t1 + t2 + t3 + t4);
        end
        tmp = [tmp(:,2), tmp(:,3), f( [x(i+3); y(:,i+3)] )];
    end
end