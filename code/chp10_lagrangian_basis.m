%% Lagrangian Linear Basis Functions
%% Figure 10.1

% Define the mesh
n = 5;               % Number of intervals for the linear basis
x = linspace(0, 1, n+1); % Nodes for linear basis
x_dense = linspace(0, 1, 500); % Dense x values for smooth plotting

% Plot linear basis functions
figure;
hold on;
for i = 1:length(x)
    % Define a linear Lagrangian basis function X_i on interval [0, 1]
    % Each basis function is triangular and centered at one of the nodes in x
    X_linear = max(1 - abs(x_dense - x(i)) / (1/(n)), 0); % Piecewise linear basis
    plot(x_dense, X_linear, 'LineWidth', 2, 'DisplayName', ['Linear Basis ', num2str(i)]);
end
title('Linear Lagrangian Basis Functions');
xlabel('x');
ylabel('X(x)');
grid on;
hold off;

