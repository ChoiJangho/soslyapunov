%% Implementation of example 1 in "Help on SOS", Andrew Pacakrd et al. 2010.
clear all;
pvar x1 x2;
pvar beta gamma;
vars = [x1; x2];
% vector of monomials
vars_mono1to6 = monomials(vars, 1:6);
vars_mono01234 = monomials(vars, [0, 1, 2, 3, 4]);
vars_mono0123 = monomials(vars, [0, 1, 2, 3]);
vars_mono012 = monomials(vars, [0, 1, 2]);
vars_mono12 = monomials(vars, [1, 2]);
vars_mono01 = monomials(vars, [0, 1]);
solver_opt.solver = 'sedumi';

eps = 1e-6;
I = eye(3);

% vector field
f = [-x2;
    x1-x2+x1^2*x2];

%% Problem 1: Fixing Lyapunov Function and solve for the maximum level.
x = vars;
A = [0, -1;
    1, -1]; % Linearization of f on x=[0,0]
Q = eye(2);
P = lyap(A', Q); % We have to transpose A, that's what the function takes.
V = x' * P * x; % Lyapunov Function
dV = diff(V,x1)*f(1)+diff(V,x2)*f(2); % lie derivative of V
l = eps * x' * x;
disp("Searching for the maximum level of the given quadratic lyapunov function with P");
disp(P);

% Bisective search on gamma.
gamma_high = 10;
gamma_low = 0;
gamma_var = (gamma_high - gamma_low) / 2 + gamma_low;
while_count = 0;
feasiblity = false;
while gamma_high - gamma_low > eps || ~feasibility
    d_gamma = gamma_high - gamma_low;
    disp(gamma_var);
    % Initialization.
    prog = sosprogram(vars);
    % Decision var: s(x)
    [prog, s] = sossosvar(prog, vars_mono01, 'wscoeff');
    % S-procedure
    t = -(l + dV) + s * (V - gamma_var);
    prog = sosineq(prog, t);
    prog = sossolve(prog,solver_opt);

    if prog.solinfo.info.pinf || prog.solinfo.info.dinf
        % Infeasible
        disp("Infeasible");
        gamma_high = gamma_var;
        feasibility = false;
    else
        % Feasible
        disp("Feasible");
        gamma_low = gamma_var;
        feasibility = true;
    end
    gamma_var = (gamma_high - gamma_low) / 2 + gamma_low;
    if while_count > 100
        disp("Infeasible in the first problem.");
        return
    end
    while_count = while_count +1;
end     

fprintf("SDP solution: max level of V (gamma) %f \n", gamma_var);

%% Problem 2: Solving for the Lyapunov Function.
h = x' * x;
l1 = eps * x' * x;
l2 = eps * x' * x;


%% Sete Iterative decision variables.
% Use the feasible solution from the previous problem. Normalize the
% Lyapunov Function.
% V_sol = V;
normalized_gamma = 10;
V_sol = normalized_gamma * V / (gamma_var - 1e-4); % 1e-4 to prevent numerical error
N_iteration = 100;
prev_beta_sol = 0;
for k = 1:N_iteration
    %% Iterative phase 1 - solving for s1, s2
    % Bisective search on beta.
    beta_low = prev_beta_sol;
    beta_high = beta_low + 10;
    beta_var = (beta_high - beta_low) / 2 + beta_low;
    feasibility = false;
    dV_sol = diff(V_sol,x1)*f(1)+diff(V_sol,x2)*f(2); % lie derivative of V    
    while_count = 0;
    while beta_high - beta_low > eps || ~feasibility
        d_beta = beta_high - beta_low;
        disp("----------");
        disp("check beta: ");
        disp(beta_var);
        progS = sosprogram(vars); % Initialization.
        % Decision var: s1(x)
        [progS, s1] = sossosvar(progS, vars_mono012, 'wscoeff');
        % Decision var: s2(x)
        [progS, s2] = sossosvar(progS, vars_mono01, 'wscoeff');

        t1 = -((beta_var - h) * s1 + (V_sol - normalized_gamma)); % (21)
%         t1 = -((beta_var - h) * s1 + (V_sol - gamma_var)); % (21)

        t2 = -(l2 + dV_sol) + s2 * (V_sol - normalized_gamma); % (22)
%         t2 = -(l2 + dV_sol) + s2 * (V_sol - gamma_var); % (22)

        progS = sosineq(progS, t1);
        progS = sosineq(progS, t2);

        progS = sossolve(progS,solver_opt);

        if progS.solinfo.info.pinf || progS.solinfo.info.dinf
            % Infeasible
            disp("Infeasible");
            beta_high = beta_var;
            feasibility = false;
        else
            % Feasible
            disp("Feasible");
            beta_low = beta_var;
            feasibility = true;
        end
        beta_var = (beta_high - beta_low) / 2 + beta_low;
        if while_count > 100
            disp("Infeasible in phase 1.");
            return
        end
        while_count = while_count +1;
    end     

    beta_sol_temp = beta_var;
    s1_sol = sosgetsol(progS,s1);
    s2_sol = sosgetsol(progS,s2);

    %% Iterative phase 2 - solving for V
    progV = sosprogram(vars); % Initialization.
    % Decision variable: V
    [progV, V] = sospolyvar(progV, vars_mono1to6, 'wscoeff');
    dV = diff(V,x1)*f(1)+diff(V,x2)*f(2); % lie derivative of V    
    % Decision variable: beta
    progV = sosdecvar(progV,beta);
    % Decision variable: gamma
    progV = sosdecvar(progV,gamma);

    t1 = -((beta - h) * s1_sol + (V - gamma)); % (21)
    t2 = -(l2 + dV) + s2_sol * (V - gamma); % (22)
    progV = sosineq(progV, V-l1);
    progV = sosineq(progV, t1);
    progV = sosineq(progV, t2);
    progV = sossetobj(progV, -beta); % objective: max beta
    progV = sossolve(progV,solver_opt);

    if progV.solinfo.info.pinf || progV.solinfo.info.dinf
        % Infeasible
        disp("Infeasible in phase 2.");
        return
    end
    V_sol = sosgetsol(progV,V);
    beta_sol = double(sosgetsol(progV, beta)) - 1e-4
    gamma_sol = double(sosgetsol(progV, gamma))
    V_sol = normalized_gamma * V_sol / (gamma_sol-1e-4);
    
    if abs(beta_sol - prev_beta_sol) < eps
        break
    end
    prev_beta_sol = beta_sol;
end
disp("-----------------FINAL SOLUTION-------------------");
fprintf("SDP solution: max level of h (beta) %f \n", beta_sol);
disp("Lyapunov function V:");
disp(V_sol);
fprintf("max level of V (gamma) %f \n", gamma_sol);
