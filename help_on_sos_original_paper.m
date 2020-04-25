%% Implementation of Example 1 in
%% Tan, Weehong, and Andrew Packard. 
%% "Stability region analysis using polynomial and composite polynomial 
%% Lyapunov functions and sum-of-squares programming." 
%% IEEE Transactions on Automatic Control 53, no. 2 (2008): 565-571.
% clear all;
pvar x1 x2 beta gamma;

x = [x1; x2];
% vector field
f = [-x2;
    x1-x2+x1^2*x2];
f_max_degree = 3; % can we get this autonomously?

n_V = 3; % half of degree of V;
n_s1 = n_V - 1; % half of degree of s1i
n_s3 = 0; % half of degree of s3i
n_A = (2*n_V - 1) + f_max_degree + n_s3; % degree of (13)
n_s2 = ceil(0.5 * n_A - n_V); % half of degree of s2
n_s0 = n_s2; % half of degree of s0

% vector of monomials
z_V = monomials(x, 2:2*n_V); % monomials in V;
% sos
z_s0 = monomials(x, 0:n_s0);
z_s1 = monomials(x, 0:n_s1); % mono vector of s1;
z_s2 = monomials(x, 0:n_s2);
z_s3 = monomials(x, 0:n_s3);

eps = 1e-6;
I = eye(3);
solver_opt.solver = 'sedumi';

%% q = 2 case
R = [0.38, -0.14; -0.14, 0.28];
l1 = 0.1 * eps * x' * x;
l2 = 0.1 * eps * x' * x;
p = x'*R*x;
A = [0, -1;
    1, -1]; % Linearization of f on x=[0,0]
Q = eye(2);
P = lyap(A', Q); % We have to transpose A, that's what the function takes.
V1_sol = x' * P * x; % Initial Lyapunov Function
V2_sol = V1_sol + l1;

N_iteration = 1000;
trace_phase2_V_coeff = cell(N_iteration, 1);
trace_beta = [];
prev_beta_sol = 0;

%% Main iteration.
% Iterate over N_iteration or until beta converges.
for k = 1:N_iteration
    %% Phase 1, optimize over s
    % Maximize over beta, bisective search on beta.
    beta_low = prev_beta_sol;
    beta_high = beta_low + 5;
    feasibility = false;
    while_count = 0;
    while beta_high - beta_low > eps || ~feasibility
        beta_var = (beta_high - beta_low) / 2 + beta_low;
        fprintf("check feasibility with beta: %.6f\n", beta_var);
        progS = sosprogram(x); % Initialization.
        % Decision vars
        [progS, s012] = sossosvar(progS, z_s0, 'wscoeff');
        [progS, s021] = sossosvar(progS, z_s0, 'wscoeff');

        [progS, s11] = sossosvar(progS, z_s0, 'wscoeff');
        [progS, s12] = sossosvar(progS, z_s0, 'wscoeff');

        [progS, s21] = sossosvar(progS, z_s0, 'wscoeff');
        [progS, s22] = sossosvar(progS, z_s0, 'wscoeff');

        [progS, s31] = sossosvar(progS, z_s0, 'wscoeff');
        [progS, s32] = sossosvar(progS, z_s0, 'wscoeff');

        
        progS = sosineq(progS, -((beta_var - p) * s11 + (V1_sol-1)));
        progS = sosineq(progS, -((beta_var - p) * s12 + (V2_sol-1))); % (12)

        dV1_sol = diff(V1_sol, x1)*f(1) + diff(V1_sol, x2)*f(2);
        dV2_sol = diff(V2_sol, x1)*f(1) + diff(V2_sol, x2)*f(2);

        eq13_1 = -((1-V1_sol)*s21 + dV1_sol * s31 + l2) - s012*(V1_sol-V2_sol);
        eq13_2 = -((1-V2_sol)*s22 + dV2_sol * s32 + l2) - s021*(V2_sol-V1_sol);
        
        progS = sosineq(progS, eq13_1);
        progS = sosineq(progS, eq13_2);
        
        progS = sossolve(progS,solver_opt); % solve feasibility problem
        s012_sol = sosgetsol(progS, s012);
        s021_sol = sosgetsol(progS, s021);
        s11_sol = sosgetsol(progS, s11);
        s12_sol = sosgetsol(progS, s12);
        s21_sol = sosgetsol(progS, s21);
        s22_sol = sosgetsol(progS, s22);
        s31_sol = sosgetsol(progS, s31);
        s32_sol = sosgetsol(progS, s32);
        
        if progS.solinfo.info.pinf || progS.solinfo.info.dinf 
            % Infeasible
            disp("Infeasible");
            beta_high = beta_var;
            feasibility = false;
        elseif progS.solinfo.info.numerr == 1 
            disp("Numerical_error");
            beta_high = beta_high - eps;
            feasibility = false;
        else
            % Feasible
            disp("Feasible");
            beta_low = beta_var;
            feasibility = true;
        end
        if while_count > 100
            disp("Infeasible in phase 1 - beta");
            return
        end
        while_count = while_count +1;
    end 

    %% Phase 2, optimize over Vi and beta
    progV = sosprogram(x);
    % Decision variables
    [progV, V1] = sospolyvar(progV, z_V);
    [progV, V2] = sospolyvar(progV, z_V);
    progV = sosdecvar(progV, beta);
    
    progV = sosineq(progV, V1-l1);
    progV = sosineq(progV, V2-l1);
    
    progV = sosineq(progV, -((beta - p) * s11_sol + (V1-1)));
    progV = sosineq(progV, -((beta - p) * s12_sol + (V2-1))); % (12)

    dV1 = diff(V1, x1)*f(1) + diff(V1, x2)*f(2);
    dV2 = diff(V2, x1)*f(1) + diff(V2, x2)*f(2);

    eq13_1 = -((1-V1)*s21_sol + dV1 * s31_sol + l2) - s012_sol*(V1-V2);
    eq13_2 = -((1-V2)*s22_sol + dV2 * s32_sol + l2) - s021_sol*(V2-V1);

    progV = sosineq(progV, eq13_1);
    progV = sosineq(progV, eq13_2);
    
    progV = sossetobj(progV, -beta);
    progV = sossolve(progV, solver_opt);
    
    beta_sol = double(sosgetsol(progV, beta));
    V1_sol = sosgetsol(progV, V1);
    V2_sol = sosgetsol(progV, V2);
    if progV.solinfo.info.pinf || progV.solinfo.info.dinf 
        % Infeasible
        disp("Infeasible in phase 2.");
        return
    elseif progV.solinfo.info.numerr == 1 || progV.solinfo.info.numerr == 2
        disp("Facing numerical issue.");
        progV = sosprogram(x);
        % Decision variables
        progV = sosdecvar(progV, gamma);
        progV = sosdecvar(progV, beta);

        progV = sosineq(progV, V1_sol-l1);
        progV = sosineq(progV, V2_sol-l1);

        progV = sosineq(progV, -((beta - p) * s11_sol + (V1_sol-gamma)));
        progV = sosineq(progV, -((beta - p) * s12_sol + (V2_sol-gamma))); % (12)

        dV1_sol = diff(V1_sol, x1)*f(1) + diff(V1_sol, x2)*f(2);
        dV2_sol = diff(V2_sol, x1)*f(1) + diff(V2_sol, x2)*f(2);

        eq13_1 = -((gamma-V1_sol)*s21_sol + dV1_sol * s31_sol + l2) - s012_sol*(V1_sol-V2_sol);
        eq13_2 = -((gamma-V2_sol)*s22_sol + dV2_sol * s32_sol + l2) - s021_sol*(V2_sol-V1_sol);

        progV = sosineq(progV, eq13_1);
        progV = sosineq(progV, eq13_2);

%         progV = sossetobj(progV, -beta);
        progV = sossolve(progV, solver_opt);

        beta_sol = 0;
        gamma_sol = sosgetsol(progV, gamma);
        V1_sol = V1_sol / gamma_sol;
        V2_sol = V2_sol / gamma_sol;
        if progV.solinfo.info.pinf || progV.solinfo.info.dinf 
            % Infeasible
            disp("Infeasible in phase 2. Cannot fix numerical issue.");
            return
        elseif progV.solinfo.info.numerr == 1 || progV.solinfo.info.numerr == 2
            disp("Cannot fix numerical issue again.");
        end
    end
        
    trace_beta = [trace_beta; beta_sol];
    if abs(beta_sol - prev_beta_sol) < eps && prev_beta_sol > 0
        break
    end
    prev_beta_sol = beta_sol;
end
