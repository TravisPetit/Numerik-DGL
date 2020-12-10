clc; clear all;
%%

%u_sol = @(x_1, x_2) abs(1/2 * x_1 + x_2) * (1/2 * x_1 + x_2) + 4 * x_1^2 * x_2;
u_sol = @(x) abs(1/2 * x(1) + x(2)) * (1/2 * x(1) + x(2)) + 4 * x(1)^2 * x(2);

%fh = @(x) -5/2 * ( (0.5 * x_1 + x_2) > 0 ) + 5/2 * ( (0.5 * x_1 + x_2) < 0 ) - 8 * x_2;
fh = @(x) -5/2 * ( (0.5 * x(1) + x(2)) > 0 ) + 5/2 * ( (0.5 * x(1) + x(2)) < 0 ) - 8 * x(2);

gh = @(x) u_sol(x);

%phi = @(x,y) -(x*y + x^2 + y^2);
%phi = @(x,y) -(x^2 + y^2);
[phi, c, xl, yl] = domains('ellipsering');
%c = -0.5;

p = 11;
%%

ns = [10,50,100,500,1000,5000];
errors = 1:6;
time_create_grid = 1:6;
time_system_matrix = 1:6;
time_system_vector = 1:6;
time_bslash_solver = 1:6;
time_count_points = 1:6;

for i = 1:6
    n = ns(i);
    xs = linspace(-1,1,n);
    ys = linspace(-1,1,n);
    
    %a)
    tic
    [C, H, G, B] = makeGrid(phi, c, xs, ys, p);
    time_create_grid(i) = toc;

    %b)
    tic
    [A, Ab] = shortleyWeller(C, H);
    time_system_matrix(i) = toc;

    %c)
    tic
    fv = evaluateOnGridDomain(fh, G);
    gv = evaluateOnGridBoundary(gh, B);
    time_system_vector(i) = toc;

    %d)
    Ab_times_gv = Ab * gv;
    tic
    b = fv - Ab_times_gv;
    time_bslash_solver(i) = toc;
    
    tic
    max(C, [], 'all');
    toc = time_count_points(i);
    
    [~, nd] = size(G);
    SOL = zeros(nd, 1);
    for q = 1:nd
        SOL(q) = u_sol(G(:,q));
    end
    
    [~, nb] = size(B);
    v = zeros(nb,1);
    for q = 1:nb
        v(q) = gh(B(:,q));
    end
    
    u = A\b;
    errors(i) = max(abs(SOL - u));
end

subplot(2,2,1)
[F, P] = makeMesh(C, H, G, B);
plotMeshFunction(F, P, u, v);
title("n = 5000, p = 11");


subplot(2,2,2)
p = loglog(ns, errors, ...
    ns, 1./ns, ...
    ns, 1./(ns.^2));
title("Convergence");
xlabel("n");
ylabel("max norm");
legend("error", "h", "h^2");
p(1).Marker = "s";
p(2).LineStyle = "--";
p(3).LineStyle = "--";

subplot(2,2,3)
p = loglog(ns, time_create_grid, ...
    ns, time_system_matrix, ...
    ns, time_system_vector, ...
    ns, time_bslash_solver, ...
    ns, time_count_points, ...
    ns, ns);
title("Timings")
xlabel("n");
ylabel("time, s");
legend("makeGrid", "system Matrix", "system Vectors", "backslash-solver", "inner points", "n");
p(1).Marker = "s";
p(2).Marker = "s";
p(3).Marker = "s";
p(4).Marker = "s";
p(5).Marker = "s";
p(6).LineStyle = "--";

subplot(2,2,4)
time = time_create_grid + time_system_matrix + time_system_vector + time_bslash_solver + time_count_points;
p = loglog(time, errors);
xlabel("time");
ylabel("error");
p.Marker = "s";
title("Work vs error");