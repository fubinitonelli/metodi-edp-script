% Eulero Indietro per Equazione del Calore

function [xh,th,uh] = chp8_backward_euler(L,N,T,K,f,uL,u0)
    h = 2*L/N; tau = T/K;
    xh = linspace(-L,L,N+1)';
    th = linspace(0,T,K+1)';

    uh = zeros(N+1,K+1);
    uh(:,1) = u0(xh);
    uh(1,:) = uL(th); uh(end,:) = uL(th);

    e = ones(2*N,1);
    A = tau/h^2*spdiags([e,-2*e,e],[-1,0,1],N-1,N-1);
    I = speye(N-1,N-1);

    for i=1:K
        F = f(xh(2:end-1),th(i+1));
        F(1) = F(1) + uL(th(i+1))/h^2; % Boundary correction
        F(end) = F(end) + uL(th(i+1))/h^2; % Boundary correction
        F = tau*F;
        uh(2:end-1,i+1) = (I-A)\(uh(2:end-1,i)+F);
    end
end