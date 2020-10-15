clc; clear all;

a = 0;
b = 12;
ya = [0; 
      0.999999;
      0.000001;
      0];

r_R_0  = 2;
r_S_0  = 0.3;
r_EE_0 = 4;
r_EI_0 = 6;
r_I_0  = 3;

r_R_1  = @(x) 2;
r_S_1  = @(x) 0.3;
r_EE_1 = @(x) 4;
r_EI_1 = @(x) 6 * (x<2) + 0.5 * (x>=2);
r_I_1  = @(x) 3;

r_R_2  = @(x) 2;
r_S_2  = @(x) 0.3;
r_EE_2 = @(x) 4 * (x<3) + 2 * (x>=3);
r_EI_2 = @(x) 6 * (x<2) + 0.5 * (x>=2);
r_I_2  = @(x) 3;

% R = z(2);
% S = z(3);
% E = z(4);
% I = z(5);
f_0 = @(z) [r_R_0 * z(5) - r_S_0 * z(2);
         r_S_0 * z(2) - r_EE_0 * z(3) * z(4) - r_EI_0 * z(3) * z(5);
         r_EE_0 * z(3) * z(4) + r_EI_0 * z(3) * z(5) - r_I_0 * z(4);
         r_I_0 * z(4) - r_R_0 * z(5)];

f_0_prime = @(z) [-r_S_0, 0, 0, r_R_0;
                  r_S_0, -r_EE_0 * z(4) - r_EI_0 * z(5), -r_EE_0 * z(3), -r_EI_0 * z(3);
                  0, r_EE_0 * z(4) + r_EI_0 * z(5), r_EE_0 * z(3) - r_I_0, r_EI_0 * z(3);
                  0, 0, r_I_0, -r_R_0];
     
f_1 = @(z) [r_R_1(z(1)) * z(5) - r_S_1(z(1)) * z(2);
         r_S_1(z(1)) * z(2) - r_EE_1(z(1)) * z(3) * z(4) - r_EI_1(z(1)) * z(3) * z(5);
         r_EE_1(z(1)) * z(3) * z(4) + r_EI_1(z(1)) * z(3) * z(5) - r_I_1(z(1)) * z(4);
         r_I_1(z(1)) * z(4) - r_R_1(z(1)) * z(5)];
     
f_1_prime = @(z) [-r_S_1(z(1)), 0, 0, r_R_1(z(1));
                  r_S_1(z(1)), -r_EE_1(z(1)) * z(4) - r_EI_1(z(1)) * z(5), -r_EE_1(z(1)) * z(3), -r_EI_1(z(1)) * z(3);
                  0, r_EE_1(z(1)) * z(4) + r_EI_1(z(1)) * z(5), r_EE_1(z(1)) * z(3) - r_I_1(z(1)), r_EI_1(z(1)) * z(3);
                  0, 0, r_I_1(z(1)), -r_R_1(z(1))];
  
f_2 = @(z) [r_R_2(z(1)) * z(5) - r_S_2(z(1)) * z(2);
         r_S_2(z(1)) * z(2) - r_EE_2(z(1)) * z(3) * z(4) - r_EI_2(z(1)) * z(3) * z(5);
         r_EE_2(z(1)) * z(3) * z(4) + r_EI_2(z(1)) * z(3) * z(5) - r_I_2(z(1)) * z(4);
         r_I_2(z(1)) * z(4) - r_R_2(z(1)) * z(5)];
     
f_2_prime = @(z) [-r_S_2(z(1)), 0, 0, r_R_2(z(1));
                  r_S_2(z(1)), -r_EE_2(z(1)) * z(4) - r_EI_2(z(1)) * z(5), -r_EE_2(z(1)) * z(3), -r_EI_2(z(1)) * z(3);
                  0, r_EE_2(z(1)) * z(4) + r_EI_2(z(1)) * z(5), r_EE_2(z(1)) * z(3) - r_I_2(z(1)), r_EI_2(z(1)) * z(3);
                  0, 0, r_I_2(z(1)), -r_R_2(z(1))];
  

% theta
theta = 0.6;
n = 2^8;
sigma = 0.5;
tol = 1e-12;
maxiter = 42;
[x_0, y_0] = thetaSchema(f_0, f_0_prime, a, b, ya, n, theta, sigma, tol, maxiter);
[x_1, y_1] = thetaSchema(f_1, f_1_prime, a, b, ya, n, theta, sigma, tol, maxiter);
[x_2, y_2] = thetaSchema(f_2, f_2_prime, a, b, ya, n, theta, sigma, tol, maxiter);

theta_0 = figure;
area(x_0.', y_0.')
title("theta 0");

theta_1 = figure;
area(x_1.', y_1.')
title("theta 1");

theta_2 = figure;
area(x_2.', y_2.')
title("theta 2");

 
%RK

alpha = [0, 1/2, 1/2, 1];
beta = [0,   0,   0, 0;
        1/2, 0,   0, 0;
        0,   1/2, 0, 0;
        0,   0,   1, 0];
gamma = [1/6, 1/3, 1/3, 1/6];
n = 2^8;

[x_0, y_0] = explizitRK(f_0, a, b, ya, n, alpha, beta, gamma);
[x_1, y_1] = explizitRK(f_1, a, b, ya, n, alpha, beta, gamma);
[x_2, y_2] = explizitRK(f_2, a, b, ya, n, alpha, beta, gamma);

rk_0 = figure;
area(x_0.', y_0.')
title("rk 0");

rk_1 = figure;
area(x_1.', y_1.')
title("rk 1");

rk_2 = figure;
area(x_2.', y_2.')
title("rk 2");


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

[x_0, y_0] = explizitRKadapt(f_0, a, b, ya, h0, alpha_DOPRI, beta_DOPRI, gamma_DOPRI, gamma_hat_DOPRI, p_DOPRI, tau, 10^(-6));
[x_1, y_1] = explizitRKadapt(f_1, a, b, ya, h0, alpha_DOPRI, beta_DOPRI, gamma_DOPRI, gamma_hat_DOPRI, p_DOPRI, tau, 10^(-6));
[x_2, y_2] = explizitRKadapt(f_2, a, b, ya, h0, alpha_DOPRI, beta_DOPRI, gamma_DOPRI, gamma_hat_DOPRI, p_DOPRI, tau, 10^(-6));

rk_a_0 = figure;
area(x_0.', y_0.')
title("rk adapt 0");

rk_a_1 = figure;
area(x_1.', y_1.')
title("rk adapt 1");

rk_a_2 = figure;
area(x_2.', y_2.')
title("rk adapt 2");
