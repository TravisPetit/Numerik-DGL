function [phi, c, xl, yl] = domains(name)
    switch name
        case 'square'
            phi = @(x, y) (x + .5) * (.5 - x) * (y + .5) * (.5 - y);
            c = 0;
            xl = [-.5001, .5001];
            yl = [-.5001, .5001];
        case 'circle'
            phi = @(x, y) - x^2 - y^2;
            c = -.5^2;
            xl = [-.5001, .5001];
            yl = [-.5001, .5001];
        case 'ellipse'
            phi = @(x, y) .5^2 - (.5*(x+y)^2 + 1.3^2*.5*(x-y)^2);
            c = 0;
            xl = [-.447, .447];
            yl = [-.447, .447];
        case 'ellipsering'
            phi = @(x, y) (.5^2 - x^2 - y^2) * (x^2 + 4*y^2 - .35^2);
            c = 0;
            xl = [-.5001, .5001];
            yl = [-.5001, .5001];
        case 'triangle'
            phi = @(x, y) (x > -.5) * (y < .5) * (x-y < 0);
            c = 0;
            xl = [-.5001, .5001];
            yl = [-.5001, .5001];
        otherwise
            rf = .5; xf = .0; yf = .0;
            rm = .3; xm = .0; ym = -.3;
            rl = .1; xl = -.2; yl = .2;
            rr = .115; xr = .2; yr = .2;
            phi = @(x, y) max((rf^2 - (x - xf)^2 - (y - yf)^2) ...
                           * ((x - xm)^2 + 3^2*(y - ym - 1.5*x^2)^2 - rm^2) ...
                           * ((x - xl)^2 + (y - yl)^2 - rl^2) ...
                           * ((x - xr)^2 + (y - yr)^2 - rr^2) ...
                           * -min(y+.1, min(.1-3*x-y, .1+3*x-y)), -1e-9);
            c = 0;
            xl = [-.5001, .5001];
            yl = [-.5001, .5001];
    end
end