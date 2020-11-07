clc; clear all;
format longG;

lambda = -5;
a = 1;
b = 3;
ya = 1;

f = @(z) lambda * z(2);
fy = @(z) lambda;

actual = @(x) exp(5 - 5 * x);
dist = @(x,y) abs(x - y);

ns = 2.^ (3:13);
ns_len = 11;

sigma = 0.5;
tol = 1e-12;
maxiter = 42;

errors_odeAB2 = zeros(ns_len,1);
errors_odeAB3 = zeros(ns_len,1);

errors_odeABM2 = zeros(ns_len,1);
errors_odeABM3 = zeros(ns_len,1);

errors_odeAM2 = zeros(ns_len,1);
errors_odeAM3 = zeros(ns_len,1);

errors_odeBDF2 = zeros(ns_len,1);
errors_odeBDF3 = zeros(ns_len,1);

for i = 1:ns_len
    n = ns(i);
    
    [~,y] = odeAB2(f, a, b, ya, n);
    errors_odeAB2(i) = dist(y(end), actual(b));
    
    [~,y] = odeAB3(f, a, b, ya, n);
    errors_odeAB3(i) = dist(y(end), actual(b));
    
    [~,y] = odeABM2(f, a, b, ya, n, 2);
    errors_odeABM2(i) = dist(y(end), actual(b));
    
    [~,y] = odeABM3(f, a, b, ya, n, 2);
    errors_odeABM3(i) = dist(y(end), actual(b));
    
    [~,y] = odeAM2(f, fy, a, b, ya, n, sigma, tol, maxiter);
    errors_odeAM2(i) = dist(y(end), actual(b));
    
    [~,y] = odeAM3(f, fy, a, b, ya, n, sigma, tol, maxiter);
    errors_odeAM3(i) = dist(y(end), actual(b));
    
    [~,y] = odeBDF2(f, fy, a, b, ya, n, sigma, tol, maxiter);
    errors_odeBDF2(i) = dist(y(end), actual(b));
    
    [~,y] = odeBDF3(f, fy, a, b, ya, n, sigma, tol, maxiter);
    errors_odeBDF3(i) = dist(y(end), actual(b));
end


subplot(2,2,1)
p = loglog(ns, errors_odeAB2, ...
    ns, errors_odeAB3, ...
    ns, 1./(10*ns.^2), ...
    ns, 1./(10*ns.^3), ...
    ns, 1./(10*ns.^4));

p(1).Marker = "s";
p(2).Marker = "s";
p(3).LineStyle = "-.";
p(4).LineStyle = "-.";
p(5).LineStyle = "-.";

title("Adams-Bashforth");
xlabel("Number of steps");
ylabel("Error at x = 3");
ylim([1e-20, 1]);
legend("AB2", "AB3", "n^{-2}","n^{-3}", "n^{-4}", 'Location','southwest');


%g = figure;
subplot(2,2,4)
p = loglog(ns, errors_odeABM2, ...
    ns, errors_odeABM3, ...
    ns, 1./(10*ns.^2), ...
    ns, 1./(10*ns.^3), ...
    ns, 1./(ns.^4));

p(1).Marker = "s";
p(2).Marker = "s";
p(3).LineStyle = "-.";
p(4).LineStyle = "-.";
p(5).LineStyle = "-.";

title("Adams-Bashforth-Moulton, m=2");
xlabel("Number of steps");
ylabel("Error at x = 3");
ylim([1e-20, 1]);
legend("ABM2", "ABM3", "n^{-2}","n^{-3}", "n^{-4}", 'Location','southwest');

%h = figure;
subplot(2,2,2)
p = loglog(ns, errors_odeAM2, ...
    ns, errors_odeAM3, ...
    ns, 1./(10*ns.^2), ...
    ns, 1./(10*ns.^3), ...
    ns, 1./(10*ns.^4));
title("Adams-Moulton");
p(1).Marker = "s";
p(2).Marker = "s";
p(3).LineStyle = "-.";
p(4).LineStyle = "-.";
p(5).LineStyle = "-.";
xlabel("Number of steps");
ylabel("Error at x = 3");
ylim([1e-20, 1]);
legend("AM2", "AM3", "n^{-2}","n^{-3}", "n^{-4}", 'Location','southwest');

subplot(2,2,3)
p = loglog(ns, errors_odeBDF2, ...
    ns, errors_odeBDF3, ...
    ns, 1./(10*ns.^2), ...
    ns, 1./(10*ns.^3), ...
    ns, 1./(10*ns.^4));
title("BDF");
p(1).Marker = "s";
p(2).Marker = "s";
p(3).LineStyle = "-.";
p(4).LineStyle = "-.";
p(5).LineStyle = "-.";
xlabel("Number of steps");
ylabel("Error at x = 3");
ylim([1e-20, 1]);
legend("BDF2", "BDF3", "n^{-2}","n^{-3}", "n^{-4}", 'Location','southwest');

sgtitle('Multistep test example, \lambda = -5')