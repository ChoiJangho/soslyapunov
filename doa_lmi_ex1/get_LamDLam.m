function [Lam, dLam] = get_LamDLam(lambda)
Lam = zeros(6);
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

dLam = zeros(6, 6, 6);
dLam(1, 6, 1) = 0.5; dLam(6, 1, 1) =0.5;
dLam(3, 3, 1) = -1;

dLam(1, 4, 2) = 0.5; dLam(4, 1, 2) = 0.5;
dLam(2, 2, 2) = -1;

dLam(1, 5, 3) = 1; dLam(5, 1, 3) = 1;
dLam(2, 3, 3) = -1; dLam(3, 2, 3) = -1;

dLam(2, 5, 4) = 1; dLam(5, 2, 4) = 1;
dLam(3, 4, 4) = -1; dLam(3, 4, 4) = -1;

dLam(2, 6, 5) = -1; dLam(6, 2, 5) = -1;
dLam(3, 5, 5) = 1; dLam(5, 3, 5) = 1;

dLam(4, 6, 6) = 0.5; dLam(6, 4, 6) = 0.5;
dLam(5, 5, 6) = -1;
end