%% Solutions to heat equation
%% Figure 2.1

f = @(x) abs(x)<=1;
p = @(x,t) 1./sqrt(4*pi*t).*exp(-(x.^2)./(4*t)); % Poisson Kernel

xh = linspace(-2*pi,2*pi,200);

plot(xh,f(xh),LineWidth=1.3);
hold on; grid on;

for i=[0.1,0.5,1]
    disp(i);
    uh = conv(f(xh),p(xh,i),"same")*(xh(2)-xh(1));
    plot(xh,uh,LineWidth=1.3);
end

legend(["f(x)","u(x,0.1)","u(x,0.5)","u(x,1)"]);
xlim([-6,6]);
ylim([-0.05,1.05]);