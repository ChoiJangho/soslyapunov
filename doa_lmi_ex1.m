%% Implementation of example 1 in "Estimation of the domain of attraction for
%% polynomial systems via LMI's", B. Tibken, CDC 2000

clear all;
pvar x1 x2 c;
vars = [x1; x2];
% vector of monomials
vars_mono2 = monomials(vars, [0, 1, 2]);
vars_mono1 = monomials(vars, [0, 1]);

eps = 1e-6;
I = eye(3);

% vector field
f = [-x1;
    -x2+x1^2*x2];

% Lyapunov Function
V = x1^2 + x2^2;
dV = diff(V,x1)*f(1)+diff(V,x2)*f(2);
r_star = 1;
c_star = 1;
c0 = c_star;

%% Step 1: check feasibility at c_star;
% prog_feas = sosprogram(vars);
% 
% % [prog, G1] = sospolymatrixvar(prog, monomials(vars, 0), [3, 3], 'symmetric');
% % [prog, G2] = sospolymatrixvar(prog, monomials(vars, 0), [3, 3], 'symmetric');
% % q1 = vars_mono1' * G1 * vars_mono1;
% % q2 = vars_mono1' * G2 * vars_mono1
% 
% % Decision var 1: q1(x)
% [prog_feas, q1] = sossosvar(prog_feas, vars_mono1, 'wscoeff');
% % Decision var 2: q2(x)
% [prog_feas, q2] = sossosvar(prog_feas, vars_mono1, 'wscoeff');
% % q0(x) at c_star, determined by formula (19)
% q0 = -dV - q1*(V-c0) - q2*(c_star-V);
% prog_feas = sosineq(prog_feas, q0-eps);
% prog_feas = sosineq(prog_feas, q1-eps);
% prog_feas = sosineq(prog_feas, q2-eps);
% solver_opt.solver = 'sedumi';
% prog_feas = sossolve(prog_feas,solver_opt);
% if prog_feas.solinfo.info.pinf || prog_feas.solinfo.info.dinf
%     disp("Infeasible solution at c_star, increase monomial degree.")
%     return
% else
%     disp("Feasible solution at c_star, proceeding to c_opt problem.")
% end
% feasible_q1 = sosgetsol(prog_feas,q1);
% feasible_q2 = sosgetsol(prog_feas,q2);


%% Step 2: solve min(eta)
eta_low = 0;
eta_high = eta_low + 100;
feasibility = false;
while_count = 0;
solver_opt.solver = 'sedumi';
 
while eta_high - eta_low > eps || ~ feasibility
    eta_var = (eta_high - eta_low) / 2 + eta_low;
    disp("check eta:"); disp(eta_var);
    % Initialization
    prog = sosprogram(vars);
    % Decision var 1: q1(x)
    [prog, q1] = sospolyvar(prog, vars_mono2, 'wscoeff');
    % Decision var 2: q2(x)
    [prog, q2] = sospolyvar(prog, vars_mono2, 'wscoeff');
    % eta*q0(x), determined by formula (19), (30)
    q0_star = -dV - q1*(V-c0) - q2*(c_star-V);
    q00 = -dV - q1*(V-c0) + q2 * V; % q11 is eta dependent terms in q0;
    q01 = -q2;

    prog = sosineq(prog, eta_var * q0_star + q01 -eps);
    prog = sosineq(prog, q0_star-eps);
    prog = sosineq(prog, q1-eps);
    prog = sosineq(prog, q2-eps);
    prog = sossolve(prog,solver_opt);
    if prog.solinfo.info.pinf || prog.solinfo.info.dinf 
        % Infeasible
        disp("Infeasible");
        eta_low = eta_var;
        feasibility = false;
    elseif prog.solinfo.info.numerr == 1 
        disp("Numerical_error");
        eta_low = eta_low + eps;
        feasibility = false;
    else
        % Feasible
        disp("Feasible");
        eta_high = eta_var;
        feasibility = true;
    end
    if while_count > 100
        disp("Infeasible in phase 1 - beta");
        return
    end
    while_count = while_count +1;
    SOL_q1 = sosgetsol(prog,q1);
    SOL_q2 = sosgetsol(prog,q2);

end
c_opt = c_star+1/eta_var;

