function [xh,uh] = chp7_centered_solver(L,N,ua,ub,f)
h = L/N;
xh = linspace(0,L,N+1)';

e = ones(N-1,1);
A = spdiags([-e 2*e -e],[-1 0 1],N-1,N-1);

F = f(xh(2:end-1));
F = (h^2)*F;
F(1) = F(1)+ua; F(end) = F(end)+ub;

uh = A\F;
uh = [ua;uh;ub];
end