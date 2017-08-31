% Compares monte carlo metropolis ising model to analytical solution
% From sweep1.fig it's clear that the critical dimensionless temperature
% lies somewhere between t = 2 and t = 2.5. Somewhat arbitrarily I take the
% critical temperature to lie exactly between these two points, t = 2.25.
% This falls roughly 2% short of the critical temperature given by the
% analytical solution t = 2/asinh(1) = 2.2692.

clc;clear;
% Run mxm ising model over temperature range t
m = 100;
t = [1:0.3:1.6 1.7:0.05:2.3 2.5:0.5:5];

% Get analytical solutions
% energy
Uonsager = @(t) -coth(2./t).*(1 + (2/pi)*(2*(tanh(2./t).^2)-1).*...
            ellipke((4*sinh(2./t).^2)./cosh(2./t).^4));

% order
Lonsager = @(t) (1 - 1./sinh(2./t).^4).^0.125;
Lo = [];
for i = 1:length(t)
    if t(i) < 2/asinh(1)
        Lo(i) = Lonsager(t(i));
    else
        Lo(i) = 0;
    end
end

% simulate
U = [];
L = [];

for i = 1:length(t)
    tic
    display(sprintf('iteration: %d of %d',i,length(t)))
    [U(i), L(i)] = ising_model_mc_met(t(i),m);
    toc
    display(sprintf('energy: %f', U(end)))
    display(sprintf('order: %f\n',L(end)))    
end

figure(2)
subplot(1,2,1)
plot(t,U,t,Uonsager(t))
title('Energy')
xlabel('t')
ylabel('normalized U')
legend({'simulated', 'analytic'})

subplot(1,2,2)
plot(t,L,t,Lo)
title('Order')
xlabel('t')
ylabel('L')
legend({'simulated', 'analytic'})