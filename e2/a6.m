% clc; clear all;
% format longG;
% 
% 
% a = 0;
% b = 10;
% ya = [0;
%       2]; 
% eta_inv = 1/0.03;
% 
% f = @(z) [eta_inv * (z(3) - z(2)^3/3 + z(2));
%          -z(2)];
% 
% fy = @(z) [eta_inv * (-z(2)^2 + 1), eta_inv;
%          -1, 0];
% fy = @(z) [eta_inv * (-z(2)^2 + 1), -1;
%           eta_inv, 0];
% 
% actual = [1.556042079300293;
%         -0.332027354872852];
%     
% dist = @(x) norm(x - actual);
% 
% ns = 2.^ (3:13);
% ns_len = 11;
% 
% sigma = 0.5;
% tol = 1e-12;
% maxiter = 42;
% 
% errors_odeAB2 = zeros(ns_len,1);
% errors_odeAB3 = zeros(ns_len,1);
% 
% errors_odeABM2 = zeros(ns_len,1);
% errors_odeABM3 = zeros(ns_len,1);
% 
% errors_odeAM2 = zeros(ns_len,1);
% errors_odeAM3 = zeros(ns_len,1);
% 
% errors_odeBDF2 = zeros(ns_len,1);
% errors_odeBDF3 = zeros(ns_len,1);
% 
% for i = 1:ns_len
%     n = ns(i);
%     
%     [~,y] = odeAB2(f, a, b, ya, n);
%     errors_odeAB2(i) = dist(y(:,end));
%     
%     [~,y] = odeAB3(f, a, b, ya, n);
%     errors_odeAB3(i) = dist(y(:,end));
%     
%     [~,y] = odeABM2(f, a, b, ya, n, 2);
%     errors_odeABM2(i) = dist(y(:,end));
%     
%     [~,y] = odeABM3(f, a, b, ya, n, 2);
%     errors_odeABM3(i) = dist(y(:,end));
%     
%     [~,y] = odeAM2(f, fy, a, b, ya, n, sigma, tol, maxiter);
%     errors_odeAM2(i) = dist(y(:,end));
%     
%     [~,y] = odeAM3(f, fy, a, b, ya, n, sigma, tol, maxiter);
%     errors_odeAM3(i) = dist(y(:,end));
%     
%     [~,y] = odeBDF2(f, fy, a, b, ya, n, sigma, tol, maxiter);
%     errors_odeBDF2(i) = dist(y(:,end));
%     
%     [~,y] = odeBDF3(f, fy, a, b, ya, n, sigma, tol, maxiter);
%     errors_odeBDF3(i) = dist(y(:,end));
% end
% 

subplot(2,2,1)
p = loglog(ns, errors_odeAB2, ...
    ns, errors_odeAB3, ...
    ns, 111600*1./(ns.^2), ...
    ns, 111600*1./(ns.^3), ...
    ns, 111600*1./(ns.^4));

p(1).Marker = "s";
p(2).Marker = "s";
p(3).LineStyle = "-.";
p(4).LineStyle = "-.";
p(5).LineStyle = "-.";

title("Adams-Bashforth");
xlabel("Number of steps");
ylabel("Error at x = 10");
ylim([1e-10, 10]);
legend("AB2", "AB3", "n^{-2}","n^{-3}", "n^{-4}", 'Location','southwest');


%g = figure;
subplot(2,2,4)
p = loglog(ns, errors_odeABM2, ...
    ns, errors_odeABM3, ...
    ns, 100000*1./(ns.^2), ...
    ns, 10000000*1./(ns.^3), ...
    ns, 1000000000*1./(ns.^4));

p(1).Marker = "s";
p(2).Marker = "s";
p(3).LineStyle = "-.";
p(4).LineStyle = "-.";
p(5).LineStyle = "-.";

title("Adams-Bashforth-Moulton, m=2");
xlabel("Number of steps");
ylabel("Error at x = 10");
ylim([1e-10, 10]);
legend("ABM2", "ABM3", "n^{-2}","n^{-3}", "n^{-4}", 'Location','southwest');

%h = figure;
subplot(2,2,2)
p = loglog(ns, errors_odeAM2, ...
    ns, errors_odeAM3, ...
    ns, 10000*1./(1*ns.^2), ...
    ns, 100000*1./(1*ns.^3), ...
    ns, 100000*1./(1*ns.^4));
title("Adams-Moulton");
p(1).Marker = "s";
p(2).Marker = "s";
p(3).LineStyle = "-.";
p(4).LineStyle = "-.";
p(5).LineStyle = "-.";
xlabel("Number of steps");
ylabel("Error at x = 10");
ylim([1e-10, 10^10]);
legend("AM2", "AM3", "n^{-2}","n^{-3}", "n^{-4}", 'Location','southwest');

subplot(2,2,3)
p = loglog(ns, errors_odeBDF2, ...
    ns, errors_odeBDF3, ...
    ns, 100000*1./(ns.^2), ...
    ns, 100000*1./(ns.^3), ...
    ns, 100000*1./(ns.^4));
title("BDF");
p(1).Marker = "s";
p(2).Marker = "s";
p(3).LineStyle = "-.";
p(4).LineStyle = "-.";
p(5).LineStyle = "-.";
xlabel("Number of steps");
ylabel("Error at x = 10");
ylim([1e-10, 10]);
legend("BDF2", "BDF3", "n^{-2}","n^{-3}", "n^{-4}", 'Location','southwest');

sgtitle('Multistep methods, Van-der-Pol: \epsilon = 0.03')