%% Implementation of example 1 in "Help on SOS", Andrew Pacakrd et al. 2010.
clear all;
pvar x1 x2;
pvar beta gamma;
vars = [x1; x2];
% vector of monomials
vars_mono1to6 = monomials(vars, 1:6);
vars_mono01234 = monomials(vars, [0, 1, 2, 3, 4]);
vars_mono123 = monomials(vars, [1, 2, 3]);
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
% Q = [5 0; 0 2];
P = lyap(A', Q); % We have to transpose A, that's what the function takes.
V = x' * P * x; % Lyapunov Function
dV = diff(V,x1)*f(1)+diff(V,x2)*f(2); % lie derivative of V
l = eps * x' * x;
disp("Searching for the maximum level of the given quadratic lyapunov function with P");
disp(P);

% Bisective search on gamma.
gamma_high = 10;
gamma_low = 0;
while_count = 0;
feasiblity = false;
while gamma_high - gamma_low > eps || ~feasibility
    gamma_var = (gamma_high - gamma_low) / 2 + gamma_low;
    disp(gamma_var);
    % Initialization.
    prog = sosprogram(vars);
    % Decision var: s(x)
    [prog, s] = sossosvar(prog, vars_mono01, 'wscoeff');
    % S-procedure
    t = -(l + dV) + s * (V - gamma_var);
    prog = sosineq(prog, t);
    prog = sossolve(prog,solver_opt);

    if prog.solinfo.info.pinf || prog.solinfo.info.dinf || prog.solinfo.info.numerr > 0
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
    if while_count > 100
        disp("Infeasible in the first problem.");
        return
    end
    while_count = while_count +1;
end     

fprintf("SDP solution: max level of V (gamma) %f \n", gamma_var);

%% Problem 2: Solving for the Lyapunov Function.
% h = x' * P * x;
h = x' * x;
l1 = 0.1 * eps * x' * x;
l2 = 0.1 * eps * x' * x;


trace_phase1_gamma = [];
trace_phase1_beta = [];
trace_phase2_gamma = [];
trace_phase2_beta = [];

%% Sete Iterative decision variables.
% Use the feasible solution from the previous problem. Normalize the
% Lyapunov Function.
% V_sol = V;
normalized_gamma = 1;
V_sol = normalized_gamma * V / (gamma_var - 1e-4); % 1e-4 to prevent numerical error
dV_sol = diff(V_sol,x1)*f(1)+diff(V_sol,x2)*f(2); % lie derivative of V    


N_iteration = 100;
trace_phase2_V_coeff = cell(100, 1);

prev_beta_sol = 0;
for k = 1:N_iteration
    %% Iterative phase 1 - solving for s1, s2
    
%     if k > 1
%     %% Maximize gamma
%     % bisective search on gamma
%     gamma_low = normalized_gamma;
%     gamma_high = gamma_low + 100;
%     gamma_var = (gamma_high - gamma_low) / 2 + beta_low;
%     feasibility = false;
%     dV_sol = diff(V_sol,x1)*f(1)+diff(V_sol,x2)*f(2); % lie derivative of V    
%     while_count = 0;  
%     while gamma_high - gamma_low > eps || ~feasibility
%         d_gamma = gamma_high - gamma_low;
%         disp("----------");
%         disp("check gamma: ");
%         disp(gamma_var);
%         progG = sosprogram(vars); % Initialization.
%         % Decision var: s1(x)
%         [progG, s1] = sossosvar(progG, vars_mono012, 'wscoeff');
%         % Decision var: s2(x)
%         [progG, s2] = sossosvar(progG, vars_mono01, 'wscoeff');
% 
%         t1 = -((prev_beta_sol - h) * s1 + (V_sol - gamma_var)); % (21)
%         t2 = -(l2 + dV_sol) + s2 * (V_sol - gamma_var); % (22)
% 
%         progG = sosineq(progG, t1);
%         progG = sosineq(progG, t2);
%         progG = sossolve(progG,solver_opt);
% 
%         if progG.solinfo.info.pinf || progG.solinfo.info.dinf
%             % Infeasible
%             disp("Infeasible");
%             gamma_high = gamma_var;
%             feasibility = false;
%         else
%             % Feasible
%             disp("Feasible");
%             gamma_low = gamma_var;
%             feasibility = true;
%         end
%         gamma_var = (gamma_high - gamma_low) / 2 + gamma_low;
%         if while_count > 100
%             disp("Infeasible in phase 1 - gamma");
%             return
%         end
%         while_count = while_count +1;
%     end
%     gamma_sol = gamma_var;
% %     s1_sol = sosgetsol(progS,s1);
% %     s2_sol = sosgetsol(progS,s2);    
%     else
%         gamma_sol = normalized_gamma;
%     end
%     
%     trace_phase1_gamma = [trace_phase1_gamma; gamma_sol];
    
    %% Maximize beta
    % Bisective search on beta.
    beta_low = prev_beta_sol;
    beta_high = beta_low + 10;
    feasibility = false;
    while_count = 0;    
    while beta_high - beta_low > eps || ~feasibility
        beta_var = (beta_high - beta_low) / 2 + beta_low;

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

    beta_sol_temp = beta_var
    s1_sol = sosgetsol(progS,s1);
    s2_sol = sosgetsol(progS,s2);

    trace_phase1_beta = [trace_phase1_beta; beta_var];

    
    %% Iterative phase 2 - solving for V
    progV = sosprogram(vars); % Initialization.
    % Decision variable: V
    [progV, V] = sospolyvar(progV, vars_mono1to6, 'wscoeff');
%     [progV, V] = sospolyvar(progV, vars_mono1to6, 'wscoeff');
    dV = diff(V,x1)*f(1)+diff(V,x2)*f(2); % lie derivative of V    
    % Decision variable: beta
    progV = sosdecvar(progV,beta);
    % Decision variable: gamma
%     progV = sosdecvar(progV,gamma);

    t1 = -((beta - h) * s1_sol + (V - normalized_gamma)); % (21)
    t2 = -(l2 + dV) + s2_sol * (V - normalized_gamma); % (22)
    progV = sosineq(progV, V-l1);
    progV = sosineq(progV, t1);
    progV = sosineq(progV, t2);
    progV = sossetobj(progV, -beta); % objective: max beta
    progV = sossolve(progV,solver_opt);

    V_sol = sosgetsol(progV,V);
    beta_sol = double(sosgetsol(progV, beta))
    dV_sol = diff(V_sol,x1)*f(1)+diff(V_sol,x2)*f(2); % lie derivative of V    

    if progV.solinfo.info.pinf || progV.solinfo.info.dinf 
        % Infeasible
        disp("Infeasible in phase 2.");
        return
    elseif progV.solinfo.info.numerr == 1
        disp("Facing numerical issue.");
        %% Define new problem to fix numerical issue.
        progV = sosprogram(vars); % Initialization.
        % Decision variable: gamma
        progV = sosdecvar(progV,gamma);
        t2 = -(l2 + dV_sol) + s2_sol * (V_sol - gamma); % (22)
        progV = sosineq(progV, t2);
        progV = sossetobj(progV, -gamma); % objective: max gamma
        progV = sossolve(progV,solver_opt);
        gamma_sol = double(sosgetsol(progV, gamma))
        V_sol = normalized_gamma * V_sol / gamma_sol;
        dV_sol = diff(V_sol,x1)*f(1)+diff(V_sol,x2)*f(2); % lie derivative of V    
        
        progV = sosprogram(vars); % Initialization.
        % Decision variable: beta
        progV = sosdecvar(progV,beta);
        t1 = -((beta - h) * s1_sol + (V_sol - normalized_gamma)); % (21)
        progV = sosineq(progV, t1);
        progV = sossetobj(progV, -beta); % objective: max beta
        progV = sossolve(progV,solver_opt);
        beta_sol = double(sosgetsol(progV, beta))        
    end
    
    trace_phase2_beta = [trace_phase2_beta; beta_sol];
%     trace_phase2_gamma = [trace_phase2_gamma; gamma_sol];
    trace_phase2_V_coeff{k} = V_sol;
    
    if abs(beta_sol - prev_beta_sol) < eps
        break
    end
    prev_beta_sol = beta_sol;
end
disp("-----------------FINAL SOLUTION-------------------");
fprintf("SDP solution: max level of h (beta) %f \n", beta_sol);
disp("Lyapunov function V:");
V_sol
fprintf("max level of V (gamma) %f \n", gamma_sol);
fprintf("max level of h (beta) %f \n", beta_sol);
