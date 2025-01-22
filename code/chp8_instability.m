%% Instabilita' del metodo Eulero Avanti, per l'equazione del calore
%% Figure 8.1

clc;
clear;

%DATA
L = pi;
T = 1;
N = 25; 
K = 25; %K = 35

% Si osserva cher per K=10 il metodo di Eulero Esplicito e' altamente
% instabile e dopo pochi passi la soluzione esplode.
ua = @(t) cos(L)*exp(t);
ub = ua;
u0 = @(x) cos(x);
f  = @(x,t) 2*cos(x).*exp(t);
u_ex = @(x,t) cos(x).*exp(t);

tau = T/K;
h   = 2*L/N;
r   = tau/(h^2);

[x,t,u_ee] = chp8_forward_euler(L,N,T,K,ua,ub,f,u0);

subplot(1,2,1)
plot(x,u_ee(:,5),LineWidth=2);
hold on; grid on;
plot(x,u_ex(x,t(5)),'--',LineWidth=2);
legend("EA","Exact")
title("t=5");

subplot(1,2,2)
plot(x,u_ee(:,20),LineWidth=2);
hold on; grid on;
plot(x,u_ex(x,t(20)),'--',LineWidth=2);
legend("EA","Exact")
title("t=20");

