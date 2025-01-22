function [p_max, p_h, nh] = chp7_error_estimate(L,N,ua1,ub,f,uex,M,scheme)
p_max = zeros(M,1); p_h = zeros(M,1);
nh = zeros(M,1);

for i=1:M
    h = L/N; nh(i) = N;
    [xh,uh] = scheme(L,N,ua1,ub,f);

    p_max(i) = max(abs(uh-uex(xh)));
    p_h(i) = sqrt(h*sum((uh-uex(xh)).^2));
    N = 2*N;
end

end
