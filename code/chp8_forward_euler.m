% Eulero Avanti per Equazione del Calore

function [xh,th,uh] = chp8_forward_euler(L,N,T,K,ua,ub,f,u0)
h = 2*L/N;
tau = T/K;

if tau > h^2/2
    fprintf("Warning! Scheme unstable\n");
end

th = linspace(0,T,K+1)';
xh = linspace(-L,L,N+1)';

uh = zeros(N+1,K+1);
uh(:,1) = u0(xh); 
uh(1,:) = ua(th); uh(end,:) = ub(th);

e = ones(2*N,1);
A = tau/h^2*spdiags([-e 2*e -e],[-1 0 1],N-1,N-1);
for k = 2:K
    F = f(xh(2:end-1),th(k));
    F(1) = F(1)+ua(th(k))/h^2;
    F(end) = F(end)+ub(th(k))/h^2;
    F = F*tau;
    uh(2:end-1,k) = (speye(N-1,N-1)-A)*uh(2:end-1,k-1)+F;
end