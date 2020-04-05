plotdomain = [-5 5 -5 5];


iteration_count = size(trace_phase2_beta, 1);


for i = 1:iteration_count
color_V = [1, 1, 1-i/iteration_count];
color_h = [1-i/iteration_count, 1, 1];
pcontour(trace_phase2_V_coeff{i}, normalized_gamma, plotdomain, 'r-'); hold on;
pcontour(h, trace_phase2_beta(i), plotdomain, 'b-');
end