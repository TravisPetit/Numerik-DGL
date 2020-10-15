clc; clear all;

digits(100);

y_lambda_eq_1 = @(x) exp(x-1);
y_lambda_eq_m_1 = @(x) exp(1-x);

f_lambda_eq_1 = @(z) z(2);
f_lambda_eq_m_1 = @(z) -z(2);

error = @(x,y) abs(x-y);

a = 1;
b = 3;
ya = 1;

exact_lambda_eq_1 = y_lambda_eq_1(b);
exact_lambda_eq_m_1 = y_lambda_eq_m_1(b);

ns = 2.^(3:13);
N = 11; %size ns

% Part a)
fy_lambda_eq_1 = @(z) 1;
fy_lambda_eq_m_1 = @(z) -1;

sigma = 0.5;
tol = 1e-12;
maxiter = 42;

thetas = [0,0.49,0.499,0.4999,0.5,0.5001,0.501,0.51,1];
T = 9; %size thetas

errors_lambda_eq_1 = zeros(N, T); % size thetas x size ns
errors_lambda_eq_m_1 = zeros(N, T);

for theta_idx = 1:T
    theta = thetas(theta_idx);
    for ns_idx = 1:N
        n = ns(ns_idx);
        
        [~, eta] = thetaSchema(f_lambda_eq_1, fy_lambda_eq_1 , a, b, ya, n, theta, sigma, tol, maxiter);
        errors_lambda_eq_1(ns_idx, theta_idx) = error(eta(end), exact_lambda_eq_1);
        
        [~, eta] = thetaSchema(f_lambda_eq_m_1, fy_lambda_eq_m_1 , a, b, ya, n, theta, sigma, tol, maxiter);
        errors_lambda_eq_m_1(ns_idx, theta_idx) = error(eta(end), exact_lambda_eq_m_1);
    end
end

part_a_lambda_eq_1 = figure;

loglog(ns, errors_lambda_eq_1, ...
    ns, 1./ns, ...
    ns, 1./(ns.^2));

legend("\theta = 0", ...
    "\theta = 0.49", ...
    "\theta = 0.499", ...
    "\theta = 0.4999", ...
    "\theta = 0.5", ...
    "\theta = 0.5001", ...
    "\theta = 0.501", ...
    "\theta = 0.51", ...
    "\theta = 1", ...
    "1/n", ...
    "1/n^2");

title("a5 a) \lambda = 1");

part_a_lambda_eq_m_1 = figure;

loglog(ns, errors_lambda_eq_m_1, ...
    ns, 1./ns, ...
    ns, 1./(ns.^2));

legend("\theta = 0", ...
    "\theta = 0.49", ...
    "\theta = 0.499", ...
    "\theta = 0.4999", ...
    "\theta = 0.5", ...
    "\theta = 0.5001", ...
    "\theta = 0.501", ...
    "\theta = 0.51", ...
    "\theta = 1", ...
    "1/n", ...
    "1/n^2");

title("a5 a) \lambda = -1");

% Part b)
alpha_euler = [0];
beta_euler = [0];
gamma_euler = [1];

alpha_heun = [0, 1];
beta_heun = [0, 0;
             1, 0];
gamma_heun = [1/2, 1/2];

alpha_classical = [0, 1/2, 1/2, 1];
beta_classical = [0,   0,   0,  0;
                  1/2, 0,   0, 0;
                  0,   1/2, 0, 0;
                  0,   0,   1, 0];
gamma_classical = [1/6, 1/3, 1/3, 1/6];

errors_lambda_eq_1_euler = zeros(N,1);
errors_lambda_eq_1_heun = zeros(N,1);
errors_lambda_eq_1_classical = zeros(N,1);

errors_lambda_eq_m_1_euler = zeros(N,1);
errors_lambda_eq_m_1_heun = zeros(N,1);
errors_lambda_eq_m_1_classical = zeros(N,1);

for ns_idx = 1:N
        n = ns(ns_idx);
        
        [~, eta] = explizitRK(f_lambda_eq_1, a, b, ya, n, alpha_euler, beta_euler, gamma_euler);
        errors_lambda_eq_1_euler(ns_idx) = error(eta(end), exact_lambda_eq_1);
        
        [~, eta] = explizitRK(f_lambda_eq_1, a, b, ya, n, alpha_heun, beta_heun, gamma_heun);
        errors_lambda_eq_1_heun(ns_idx) = error(eta(end), exact_lambda_eq_1);
        
        [~, eta] = explizitRK(f_lambda_eq_1, a, b, ya, n, alpha_classical, beta_classical, gamma_classical);
        errors_lambda_eq_1_classical(ns_idx) = error(eta(end), exact_lambda_eq_1);
        
        % ---
        
        [~, eta] = explizitRK(f_lambda_eq_m_1, a, b, ya, n, alpha_euler, beta_euler, gamma_euler);
        errors_lambda_eq_m_1_euler(ns_idx) = error(eta(end), exact_lambda_eq_m_1);
        
        [~, eta] = explizitRK(f_lambda_eq_m_1, a, b, ya, n, alpha_heun, beta_heun, gamma_heun);
        errors_lambda_eq_m_1_heun(ns_idx) = error(eta(end), exact_lambda_eq_m_1);
        
        [~, eta] = explizitRK(f_lambda_eq_m_1, a, b, ya, n, alpha_classical, beta_classical, gamma_classical);
        errors_lambda_eq_m_1_classical(ns_idx) = error(eta(end), exact_lambda_eq_m_1);
end

part_b_lambda_eq_1 = figure;

loglog(ns, errors_lambda_eq_1_euler, ...
    ns, errors_lambda_eq_1_heun, ...
    ns, errors_lambda_eq_1_classical, ...
    ns, 1./ns, ...
    ns, 1./(ns.^2), ...
    ns, 1./(ns.^4));
legend("ERK Euler", "ERK Heun", "ERK vanilla", "1/n", "1/n^2", "1/n^4");
title("a5 b) \lambda = 1");


part_b_lambda_eq_m_1 = figure;

loglog(ns, errors_lambda_eq_m_1_euler, ...
    ns, errors_lambda_eq_m_1_heun, ...
    ns, errors_lambda_eq_m_1_classical, ...
    ns, 1./ns, ...
    ns, 1./(ns.^2), ...
    ns, 1./(ns.^4));
legend("ERK Euler", "ERK Heun", "ERK vanilla", "1/n", "1/n^2", "1/n^4");
title("a5 b) \lambda = -1");

% Part c)
alpha_RKF = [0, 1, 1/2];
beta_RKF = [0,    0,  0;
            1,    0,  0;
            1/4, 1/4, 0];
gamma_RKF = [1/2, 1/2, 0];
gamma_hat_RKF = [1/6, 1/6, 4/6];

alpha_DOPRI = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];

beta_DOPRI = [
0,          0,           0,          0,        0,           0;
1/5,        0,           0,          0,        0,           0;
3/40,       9/40,        0,          0,        0,           0;
44/45,      -56/15, 	 32/9,       0,        0,           0;
19372/6561, -25360/2187, 64448/6561, -212/729, 0,           0;
9017/3168,  -355/33, 	 46732/5247, 49/176    -5103/18656, 0;
35/384,     0,           500/1113, 	 125/192,  -2187/6784,  11/84];

gamma_DOPRI = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
gamma_hat_DOPRI= [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];

p_RKF = 2;
p_DOPRI = 4;

epsilons = 2.^(-(1:20));
E = 20; % size epsilons
tau = 0.9;

errors_RKF_lambda_eq_1 = zeros(E, 1);
errors_RKF_lambda_eq_m_1 = zeros(E, 1);

errors_DOPRI_lambda_eq_1 = zeros(E, 1);
errors_DOPRI_lambda_eq_m_1 = zeros(E, 1);

h0 = 0.5;

for eps = 1:E
    epsilon = epsilons(eps);
    
    [~, eta] = explizitRKadapt(f_lambda_eq_1, a, b, ya, h0, alpha_DOPRI, beta_DOPRI, gamma_DOPRI, gamma_hat_DOPRI, p_DOPRI, tau, epsilon);
    errors_DOPRI_lambda_eq_1(eps) = error(eta(end), exact_lambda_eq_1);
    
    [~, eta] = explizitRKadapt(f_lambda_eq_1, a, b, ya, h0, alpha_RKF, beta_RKF, gamma_RKF, gamma_hat_RKF, p_RKF, tau, epsilon);
    errors_RKF_lambda_eq_1(eps) = error(eta(end), exact_lambda_eq_1);

    % ---
    
    [~, eta] = explizitRKadapt(f_lambda_eq_m_1, a, b, ya, h0, alpha_DOPRI, beta_DOPRI, gamma_DOPRI, gamma_hat_DOPRI, p_DOPRI, tau, epsilon);
    errors_DOPRI_lambda_eq_m_1(eps) = error(eta(end), exact_lambda_eq_m_1);
    
    [~, eta] = explizitRKadapt(f_lambda_eq_m_1, a, b, ya, h0, alpha_RKF, beta_RKF, gamma_RKF, gamma_hat_RKF, p_RKF, tau, epsilon);
    errors_RKF_lambda_eq_m_1(eps) = error(eta(end), exact_lambda_eq_m_1);
    
end

part_c_lambda_eq_1 = figure;

loglog(epsilons, errors_RKF_lambda_eq_1, ...
    epsilons, errors_DOPRI_lambda_eq_1, ...
    epsilons, epsilons.^(2/3), ...
    epsilons, epsilons.^(4/5));
legend("error RKF", "error DOPRI", "\epsilon^{2/3}", "\epsilon^{4/5}");
title("a5 c) \lambda = 1");

part_c_lambda_eq_m_1 = figure;
loglog(epsilons, errors_RKF_lambda_eq_m_1, ...
    epsilons, errors_DOPRI_lambda_eq_m_1, ...
    epsilons, epsilons.^(2/3), ...
    epsilons, epsilons.^(4/5));
legend("error RKF", "error DOPRI", "\epsilon^{2/3}", "\epsilon^{4/5}");
title("a5 c) \lambda = -1");