%% Mollifiers: plot of gauss and poisson kernel
%% Figure 1.2

gauss = @(x) exp(-pi*abs(x).^2);
poisson = @(x) 1/pi*1./(1+x.^2);

xh = linspace(-2,2,1000);

subplot(1,2,1);
plot(xh,gauss(xh),"Color","blue","LineWidth",1.5);
title("Gaussian Kernel");
grid on; 
subplot(1,2,2);
plot(xh,poisson(xh),"Color","red","LineWidth",1.5);
title("Poisson Kernel");
grid on;