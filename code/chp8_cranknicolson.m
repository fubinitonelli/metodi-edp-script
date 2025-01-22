% Crank-Nicolson per Equazione del Calore

function [xh,th,uh] = chp8_cranknicolson(L,N,T,K,f,uL,u0)
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
        Fa = f(xh(2:end-1),th(i));
        Fa(1) = Fa(1)+uL(th(i))/h^2; Fa(end) = Fa(end)+uL(th(i))/h^2;
        
        Fb = f(xh(2:end-1),th(i+1));
        Fb(1) = Fb(1)+uL(th(i+1))/h^2; Fb(end) = Fb(end)+uL(th(i+1))/h^2;
        
        uh(2:end-1,i+1) = (I-0.5*A)\((I+0.5*A)*uh(2:end-1,i)+0.5*tau*Fa+0.5*tau*Fb);
    end
end