% Example of how to use ODE23 in MATLAB:  make two files
% such as call_plant.m and plant.m as shown below;
% note that in this example, the differential equations 
% in plant.m contain a switch function.

% At the MATLAB prompt, type: "call_plant" <br>
% <br>
% call_plant.m <br>
N_init_test = 25;
perturbation = 4;
T_vec = [];
x_vec = [];
x0_vec = [];

for i=1:N_init_test
    
x0 = [perturbation*rand() - 0.5*perturbation, 0.25 * (perturbation*rand() - 0.5*perturbation)];    
[T,x]=ode23(@plant, [0 500], x0);
T_vec = [T_vec; T];
x_vec = [x_vec; x];
c = linspace(1,10,size(x, 1));
plot(x(:, 1), x(:, 2)); hold on;

x0_vec = [x0_vec; x0];
end
grid on;
plot(x0_vec(:, 1), x0_vec(:, 2), 'LineStyle', 'none', 'Marker', 'd', 'Color', 'r'); hold off;


% plant.m <br>
function dx = plant(t,x)
dx1 = -x(2);
dx2 = x(1) - x(2) + x(1)^2 * x(2);

dx = [dx1; dx2];

end
