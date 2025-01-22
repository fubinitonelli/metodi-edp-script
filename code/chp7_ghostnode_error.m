%% Valutazione dell'ordine di convergenza quando si risolve 
%% il problema modello discretizzato con la condizione di Neumann
%% Figure 7.1

clear;
close all;
clc;

disp('Es 2.2')
disp('Valutazione dell''ordine di convergenza quando si risolve')
disp('il problema modello discretizzato con la condizione di Neumann')

N = 100;
L = 1;
ua1 = 0;
ub  = 0;
f = @(x) 4*pi^2*cos(2*pi*x);
M = 6;

% Soluzione esatta
uex = @(x) cos(2*pi*x)-1;

% Valutazione dell'ordine di convergenza p (risp. in norma del massimo e
% norma h)

% Upwind
upwind = @(a,b,c,d,e) chp7_upwind_solver(a,b,c,d,e);
centrato = @(a,b,c,d,e) chp7_centered_solver(a,b,c,d,e);

[e_max_up, e_h_up, nh] = chp7_error_estimate(L,N,ua1,ub,f,uex,M,upwind);
[e_max_ce, e_h_ce, ~] = chp7_error_estimate(L,N,ua1,ub,f,uex,M,centrato);

figure

subplot(1,2,1);
loglog(nh,e_max_up,'--','LineWidth',2);
hold on; grid on;
loglog(nh,e_h_up,'-o','LineWidth',2);
loglog(nh,1./nh,'LineWidth',2);
legend(["E_{inf}","E_h","h"]);
title("Upwind Error");

subplot(1,2,2);
loglog(nh,e_max_ce,'--','LineWidth',2);
hold on; grid on;
loglog(nh,e_h_ce,'-o','LineWidth',2);
loglog(nh,1./nh.^2,'LineWidth',2);
legend(["E_{inf}","E_h","h^2"]);
title("Ghost Node Error");