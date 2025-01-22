% Upwind solver for mixed boundary condtion
% This is a first order method, so we expect degradation in  
% the precision of the solution. 

function [xh,uh] = chp7_upwind_solver(L,N,ua1,ub,f)
    xh = linspace(0,L,N+1)';
    h = L/N;

    % Costruzione della matrice A
    e = ones(N,1);
    A = spdiags([-e 2*e -e],[-1 0 1],N,N);
    A(1,1) = 1; A(1,2) = -1; % Correzione di Neumann

    % Costruzione del termine noto f
    F = f(xh(1:end-1));
    F = h^2*F;

    % Correzione termine noto con termine condizione al bordo
    F(1) = h*ua1; % Neumann 
    F(end) = F(end)+ub; % Dirichlet

    % Risoluzione del sistema lineare
    uh = A\F;
    uh = [uh;ub];
end