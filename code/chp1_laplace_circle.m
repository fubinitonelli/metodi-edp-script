%% Solving analytically Laplace equation on the circle
%% Figure 1.1

f = @(theta) theta.*(pi-theta);
r = 1;

% Define Circular Domain
N = 40; % number of points
theta = linspace(0, 2 * pi, N);
R = linspace(0, r, N);
[thetaG, RG] = meshgrid(theta, R);
X = RG .* cos(thetaG);
Y = RG .* sin(thetaG);
[m, n] = size(X);

% Fourier coefficients for n = 0
a0 = FourierCoeff(f, 2 * pi, 0);
Z = a0 / 2 * ones(m, n);

% Add higher-order terms
for i = 1:5
    [a, b] = FourierCoeff(f, 2 * pi, i);
    Z = Z + (RG / r) .^ i .* (a * cos(i * thetaG) + b * sin(i * thetaG));
end

% Surf plot
figure;
surf(X, Y, Z,'FaceAlpha',0.5);
colormap('parula');
%colorbar;
title(sprintf('Equazione di Laplace sul Cerchio r = %d', r));
xlabel('X');
ylabel('Y');
zlabel('Z');

% Plot boundary circle on XY plane
hold on;
circle_x = r * cos(theta);
circle_y = r * sin(theta);
plot3(circle_x, circle_y, zeros(size(circle_x)), 'b-', 'LineWidth', 2);
hold off;
