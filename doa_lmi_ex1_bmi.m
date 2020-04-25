%% Implementation of example 1 in "Estimation of the domain of attraction for
%% polynomial systems via LMI's", B. Tibken, CDC 2000
%% Solving Bilinear SDP Using BMI Solver (based on)
%% Dinh, Quoc Tran, Suat Gumussoy, Wim Michiels, and Moritz Diehl. 
%% "Combining convex?concave decompositions and linearization approaches 
%% for solving BMIs, with application to static output feedback." 
%% IEEE Transactions on Automatic Control 57, no. 6 (2011): 1377-1390.
clear all;

% constants
c_star = 1;
c0 = c_star;
eps = 1e-8;

%% Phase1
% Decision variable
G1 = sdpvar(3, 3, 'symmetric'); % q1 Gram Mtx
G2 = sdpvar(3, 3, 'symmetric'); % q2 Gram Mtx
lambda = sdpvar(6, 1);

eta0 = 10;
% V = x1^2 + x2^2
% Belows are Gram Mtx representation. z = [1 x1 x2 x1^2 x1x2 x2^2]
dV = diag([0, -2, -2, 0, 2, 0]);
q1 = double2sdpvar(zeros(6)); q1(1:3, 1:3) = G1;
q2 = double2sdpvar(zeros(6)); q2(1:3, 1:3) = G2;

q1V = double2sdpvar(zeros(6));
q1V([2,4,5],[2,4,5]) = G1;
q1V([3,5,6],[3,5,6]) = q1V([3,5,6],[3,5,6]) + G1;

q2V = double2sdpvar(zeros(6));
q2V([2,4,5],[2,4,5]) = G2;
q2V([3,5,6],[3,5,6]) = q2V([3,5,6],[3,5,6]) + G2;

Lam = double2sdpvar(zeros(6));
Lam(1, 4) = 0.5 * lambda(2); Lam(4, 1) = Lam(1, 4);
Lam(1, 5) = lambda(3); Lam(5, 1) = Lam(1, 5);
Lam(1, 6) = 0.5 * lambda(1); Lam(6, 1) = Lam(1, 6);
Lam(2, 2) = -lambda(2);
Lam(2, 3) = -lambda(3); Lam(3, 2) = Lam(2, 3);
Lam(2, 5) = lambda(4); Lam(5, 2) = Lam(2, 5);
Lam(2, 6) = -lambda(5); Lam(6, 2) = Lam(2, 6);
Lam(3, 3) = -lambda(1);
Lam(3, 4) = -lambda(4); Lam(4, 3) = Lam(3, 4);
Lam(3, 5) = lambda(5); Lam(5, 3) = Lam(3, 5);
Lam(4, 6) = 0.5 * lambda(6); Lam(6, 4) = Lam(4, 6);
Lam(5, 5) = -lambda(6);

q0_star = -dV - q1V + c0 * q1 - c_star * q2 + q2V + Lam;
epsI6 = eps * eye(6); % PSD -> PD Margin
epsI3 = eps * eye(3);
% Constraints
Fx = [-q2+eta0*q0_star-epsI6 >= 0, ...
    q0_star-epsI6 >=0, ...
    G1-epsI3 >=0, G2-epsI3 >=0];
diagnostic = optimize(Fx);
if diagnostic.problem < 0
    disp("Feasibility problem with initial eta0 failed!");
else
    G1_val = double(G1);
    G2_val = double(G2);
    lambda_val = double(lambda);
end

%% Phase2
% initialize x_k
etak = eta0;
G1k = G1_val;
G2k = G2_val;
lambdak = lambda_val;

% Hyperparameters. Can we make it adaptive?
rho = 1;
Q_G1 = 1e-2;
Q_G2 = 1e-2;
Q_lambda = 1e-2;
Q_eta = 1;

iterations = 1000;
for i = 1:iterations
    prev_etak = etak;
    G1kl = zeros(6, 1);
    G1kl(1) = G1k(1, 1); G1kl(2) = G1k(1, 2); G1kl(3) = G1k(1, 3);
    G1kl(4) = G1k(2, 2); G1kl(5) = G1k(2, 3); G1kl(6) = G1k(3, 3);

    G2kl = zeros(6, 1);
    G2kl(1) = G2k(1, 1); G2kl(2) = G2k(1, 2); G2kl(3) = G2k(1, 3);
    G2kl(4) = G2k(2, 2); G2kl(5) = G2k(2, 3); G2kl(6) = G2k(3, 3);

    
    % initialize decision variables for local opt. problem.
    G1 = sdpvar(3, 3, 'symmetric'); % q1 Gram Mtx
    G2 = sdpvar(3, 3, 'symmetric'); % q2 Gram Mtx
    lambda = sdpvar(6, 1);
    eta = sdpvar(1, 1);
    
    dV = diag([0, -2, -2, 0, 2, 0]);
    q1 = double2sdpvar(zeros(6)); q1(1:3, 1:3) = G1;
    q2 = double2sdpvar(zeros(6)); q2(1:3, 1:3) = G2;

    q1V = double2sdpvar(zeros(6));
    q1V([2,4,5],[2,4,5]) = G1;
    q1V([3,5,6],[3,5,6]) = q1V([3,5,6],[3,5,6]) + G1;
    q2V = double2sdpvar(zeros(6));
    q2V([2,4,5],[2,4,5]) = G2;
    q2V([3,5,6],[3,5,6]) = q2V([3,5,6],[3,5,6]) + G2;

    Lam = double2sdpvar(zeros(6));
    Lam(1, 4) = 0.5 * lambda(2); Lam(4, 1) = Lam(1, 4);
    Lam(1, 5) = lambda(3); Lam(5, 1) = Lam(1, 5);
    Lam(1, 6) = 0.5 * lambda(1); Lam(6, 1) = Lam(1, 6);
    Lam(2, 2) = -lambda(2);
    Lam(2, 3) = -lambda(3); Lam(3, 2) = Lam(2, 3);
    Lam(2, 5) = lambda(4); Lam(5, 2) = Lam(2, 5);
    Lam(2, 6) = -lambda(5); Lam(6, 2) = Lam(2, 6);
    Lam(3, 3) = -lambda(1);
    Lam(3, 4) = -lambda(4); Lam(4, 3) = Lam(3, 4);
    Lam(3, 5) = lambda(5); Lam(5, 3) = Lam(3, 5);
    Lam(4, 6) = 0.5 * lambda(6); Lam(6, 4) = Lam(4, 6);
    Lam(5, 5) = -lambda(6);
    q0_star = -dV - q1V + c0 * q1 - c_star * q2 + q2V + Lam;

    % vectorization
    G1l = double2sdpvar(zeros(6, 1));
    G1l(1) = G1(1, 1); G1l(2) = G1(1, 2); G1l(3) = G1(1, 3);
    G1l(4) = G1(2, 2); G1l(5) = G1(2, 3); G1l(6) = G1(3, 3);
    G2l = double2sdpvar(zeros(6, 1));
    G2l(1) = G2(1, 1); G2l(2) = G2(1, 2); G2l(3) = G2(1, 3);
    G2l(4) = G2(2, 2); G2l(5) = G2(2, 3); G2l(6) = G2(3, 3); 
    
    % Eq (4.1) in the paper
    f_lin = eta;
    f_lin = f_lin + 0.5 * rho * Q_G1 * trace((G1 - G1k)' * (G1-G1k)); 
    f_lin = f_lin + 0.5 * rho * Q_G2 * trace((G2 - G2k)' * (G2-G2k)); 
    f_lin = f_lin + 0.5 * rho * Q_lambda * trace((lambda - lambdak)' * (lambda-lambdak)); 
    f_lin = f_lin + 0.5 * rho * Q_eta * (eta - etak)^2;
    [~, H, dH_eta, dH_G1, dH_G2, dH_Lam] = linearizeConstraint(etak, G1k, G2k, lambdak);
%     G = 0.25 * (eta * eye(6) - q0_star)' * (eta * eye(6) - q0_star);

    % Constraints
    % Linearized Main Constraint.
    Hk = H + dH_eta * (eta - etak);    
    for j = 1:6
        Hk = Hk + (G1l(j) - G1kl(j)) * dH_G1(:, :, j); 
        Hk = Hk + (G2l(j) - G2kl(j)) * dH_G2(:, :, j); 
        Hk = Hk + (lambda(j) - lambdak(j)) * dH_Lam(:, :, j); 
    end
    linconstr = [-q2+Hk-epsI6, 0.5*(eta * eye(6) - q0_star);
        0.5*(eta * eye(6) - q0_star)', eye(6)];
    Fx = [linconstr >= 0, ...
        q0_star-epsI6 >=0, ...
        G1-epsI3 >=0, G2-epsI3 >=0];
    diagnostic = optimize(Fx, f_lin);
    if diagnostic.problem < 0
        disp("Iteration failed!");
    else
        etak = double(eta)
        G1kl = double(G1);
        G2kl = double(G2);
        lambdak = double(lambda);
    end
    
    if abs(prev_etak - etak) < eps
        fprintf("Iteration converged at count=%d.\n", i);
        fprintf("eta* = %.4f\n", etak);
        break
    end
end


function [G, H, dH_eta, dH_G1, dH_G2, dH_Lam] = linearizeConstraint(eta, G1, G2, lambda)
% eta, G1, G2, lambda are all double(array) (not sdpvar).
% Linearizing the concave part.
% H = (eta * I - q0_star)' * (eta * I - q0_star)
c_star = 1;
c0 = c_star;
dV = diag([0, -2, -2, 0, 2, 0]);

q1 = zeros(6); q1(1:3, 1:3) = G1;
q2 = zeros(6); q2(1:3, 1:3) = G2;

q1V = zeros(6);
q1V([2,4,5],[2,4,5]) = G1;
q1V([3,5,6],[3,5,6]) = q1V([3,5,6],[3,5,6]) + G1;

q2V = zeros(6);
q2V([2,4,5],[2,4,5]) = G2;
q2V([3,5,6],[3,5,6]) = q2V([3,5,6],[3,5,6]) + G2;

[Lam, dLam] = get_LamDLam(lambda);
q0_star = -dV - q1V + c0 * q1 - c_star * q2 + q2V + Lam;
dq = get_dq();
dqV = get_dqV();
dq0_star_G1 = -dqV + c0 * dq;
dq0_star_G2 = dqV - c_star * dq;
dq0_star_lam = dLam;
                   
G = 0.25 * (eta * eye(6) - q0_star)' * (eta * eye(6) - q0_star);
H = 0.25 * (eta * eye(6) + q0_star)' * (eta * eye(6) + q0_star);
dH_eta = 0.25 * ((eta * eye(6) + q0_star)' + (eta * eye(6) + q0_star));

dH_G1 = zeros(6,6,6);
dH_G2 = zeros(6,6,6);
dH_Lam = zeros(6,6,6);
for i = 1:6
dH_G1(:,:,i) = 0.25 * ((eta * eye(6) + q0_star)' * dq0_star_G1(:, :, i) ...
    + dq0_star_G1(:, :, i)' * (eta * eye(6) + q0_star));
dH_G2(:,:,i) = 0.25 * ((eta * eye(6) + q0_star)' * dq0_star_G2(:, :, i) ...
    + dq0_star_G2(:, :, i)' * (eta * eye(6) + q0_star));    
dH_Lam(:,:,i) = 0.25 * ((eta * eye(6) + q0_star)' * dq0_star_lam(:, :, i) ...
    + dq0_star_lam(:, :, i)' * (eta * eye(6) + q0_star));
end
end


