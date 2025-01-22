%% Plotting the analytical solution to the homogeneous advection-diffusion equation
%% Figure 11.1

xh = linspace(0,1,1000);

% epsilon fixed at one 

uex = @(x,beta) (exp(beta*x)-1)/(exp(beta)-1); % exact solution

% Plotting for different values of beta, the bigger it gets,
% the more the advection term is dominant

plot(xh,uex(xh,0.00001))
hold on; grid on;
plot(xh,uex(xh,3))
plot(xh,uex(xh,25))

legend("\beta/\epsilon = 0","\beta/\epsilon = 3","\beta/\epsilon = 25",Location="northwest");
