%% Difference Equation Solution - Oscillating Behaviour
%% Figure 11.2

clc;
clear; 
close all;

N = 50; % number of sub-intervals
xv = linspace(0,1,N+1);

beta = 13; epsilon = 0.6;
uex = @(x) (exp(beta/epsilon*x)-1)/(exp(beta/epsilon)-1); % exact solution

% Difference Equation Solution
N = 5; h = 1/N;
xh = linspace(0,1,N+1);
Pe = beta*h/(2*epsilon); % Peclet's Number

fprintf("Peclet Number = %d\n",Pe);

A_1 = 1/(((1+Pe)/(1-Pe))^N-1); 
A_2 = -A_1;
rho = (1+Pe)/(1-Pe);

uh = zeros(N+1,1);
for i=0:N
    uh(i+1) = A_1*rho^i+A_2;
end

plot(xv,uex(xv),'LineWidth',1);
hold on; grid on;
plot(xh,uh,'-o','LineWidth',1);
legend("u_{ex}(x)","Pe = 2.16",location="northwest");