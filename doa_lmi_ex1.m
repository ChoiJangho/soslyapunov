%% Implementation of example 1 in "Estimation of the domain of attraction for
%% polynomial systems via LMI's", B. Tibken, CDC 2000

clear all;
pvar x1 x2 eta c;
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
prog_feas = sosprogram(vars);

% [prog, G1] = sospolymatrixvar(prog, monomials(vars, 0), [3, 3], 'symmetric');
% [prog, G2] = sospolymatrixvar(prog, monomials(vars, 0), [3, 3], 'symmetric');
% q1 = vars_mono1' * G1 * vars_mono1;
% q2 = vars_mono1' * G2 * vars_mono1

% Decision var 1: q1(x)
[prog_feas, q1] = sossosvar(prog_feas, vars_mono1, 'wscoeff');
% Decision var 2: q2(x)
[prog_feas, q2] = sossosvar(prog_feas, vars_mono1, 'wscoeff');
% q0(x) at c_star, determined by formula (19)
q0 = -dV - q1*(V-c0) - q2*(c_star-V);
prog_feas = sosineq(prog_feas, q0-eps);
prog_feas = sosineq(prog_feas, q1-eps);
prog_feas = sosineq(prog_feas, q2-eps);
solver_opt.solver = 'sedumi';
prog_feas = sossolve(prog_feas,solver_opt);
if prog_feas.solinfo.info.pinf || prog_feas.solinfo.info.dinf
    disp("Infeasible solution at c_star, increase monomial degree.")
    return
else
    disp("Feasible solution at c_star, proceeding to c_opt problem.")
end
feasible_q1 = sosgetsol(prog_feas,q1);
feasible_q2 = sosgetsol(prog_feas,q2);


%% Step 2: solve min(eta)
% Initialization
prog = sosprogram(vars);
% 
% % Decision var 1: q0(x)
% [prog, q0] = sossosvar(prog, vars_mono1, 'wscoeff');
% % Decision var 2: q1(x)
% [prog, q1] = sossosvar(prog, vars_mono1, 'wscoeff');
% % Decision var 3: q2(x)
% [prog, q2] = sossosvar(prog, vars_mono1, 'wscoeff');
% 
% % prog = sosdecvar(prog,c);
% 
% prog = soseq(prog, q0 + q1*(V-c0) + q2*(c-V) +dV);
% 
% prog = sossetobj(prog, -c);
% 
% solver_opt.solver = 'sedumi';
% prog = sossolve(prog,solver_opt);
% SOL_q1 = sosgetsol(prog,q1);
% SOL_q2 = sosgetsol(prog,q2);

%%%%%%%%%%%%
% [prog, G1] = sospolymatrixvar(prog, monomials(vars, 0), [3, 3], 'symmetric');
% [prog, G2] = sospolymatrixvar(prog, monomials(vars, 0), [3, 3], 'symmetric');
% [prog, G00] = sospolymatrixvar(prog, monomials(vars, 0), [3, 3], 'symmetric');
% [prog, G01] = sospolymatrixvar(prog, monomials(vars, 0), [3, 3], 'symmetric');
prog = sosdecvar(prog,eta);
% 
% q1 = vars_mono1' * G1 * vars_mono1;
% q2 = vars_mono1' * G2 * vars_mono1;
% q00 = vars_mono1' * G00 * vars_mono1;
% q01= vars_mono1' * G01 * vars_mono1;
% 
% prog = soseq(prog, q00 + dV + q1*(V-c0) + q2 * c_star);
% prog = soseq(prog, q01 + q2 * (1-V));
% % I cannot apply this inequality condition.
% prog = sosmatrixineq(prog, eta*G00+G01 - eps*I, 'quadraticMineq');
% prog = sosmatrixineq(prog,G00-eps*I,'quadraticMineq');
% prog = sosmatrixineq(prog,G1-eps*I,'quadraticMineq');
% prog = sosmatrixineq(prog,G2-eps*I,'quadraticMineq');


%%%%%%%%%%%%%
% Decision var 1: q1(x)
[prog, q1] = sossosvar(prog, vars_mono1, 'wscoeff');
% Decision var 2: q2(x)
[prog, q2] = sossosvar(prog, vars_mono1, 'wscoeff');


% eta*q0(x), determined by formula (19), (30)
etaq0 = -eta * dV - eta * q1*(V-c0) - q2*(eta * c_star + 1 -V);
q00 = -dV - q1*(V-c0) - q2 * c_star; % q11 is eta dependent terms in q0;
q01 = -q2 * (1-V);
prog = sosineq(prog, etaq0-eps);
prog = sosineq(prog, q00-eps);
prog = sosineq(prog, q1-eps);
prog = sosineq(prog, q2-eps);


% set objective: min eta
prog = sossetobj(prog, eta);

solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);
SOL_q1 = sosgetsol(prog,q1);
SOL_q2 = sosgetsol(prog,q2);
