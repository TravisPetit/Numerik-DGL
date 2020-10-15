clc; clear all;

epsilon = 0.03;
eps_inv = 1/epsilon;

a = 0;
b = 10;
ya = [0;2];

f = @(z) [eps_inv * (z(3) - z(2)^3/3 + z(2));
          -z(2)];

fy = @(z) [eps_inv * (z(2)^2 + 1), eps_inv;
           -1,                         0];

% theta
theta = 0.6;
n = 2^8;
sigma = 0.5;
tol = 1e-12;
maxiter = 42;
[x, y] = thetaSchema(f, fy, a, b, ya, n, theta, sigma, tol, maxiter);
theta_f = figure;
plot(x,y)
title("theta");

% RK
alpha = [0, 1/2, 1/2, 1];
beta = [0,   0,   0, 0;
        1/2, 0,   0, 0;
        0,   1/2, 0, 0;
        0,   0,   1, 0];
gamma = [1/6, 1/3, 1/3, 1/6];
[x, y] = explizitRK(f, a, b, ya, n, alpha, beta, gamma);

rk_f = figure;
plot(x, y)
title("rk");

% RK adapt
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
p_DOPRI = 4;

h0 = 0.1;
tau = 0.9;

[x, y] = explizitRKadapt(f, a, b, ya, h0, alpha_DOPRI, beta_DOPRI, gamma_DOPRI, gamma_hat_DOPRI, p_DOPRI, tau, 10^(-6));

rk_adapt_f = figure;
plot(x, y);
title("rk adapt");
