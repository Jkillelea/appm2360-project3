clear; clc; close all;

% params
alpha = 0.00001;
beta  = 0.00003;
gamma = 0.00001;
N0    = 50000;

% SZR Model
dS = @(S, Z) -beta.*Z.*S;
dZ = @(S, Z) (beta - alpha).*Z.*S + gamma.*(N0 - S - Z);

[S, Z] = meshgrid(0:(N0/100):N0, 0:(N0/100):N0);

data_ds = dS(S, Z);
data_dz = dZ(S, Z);

fig = figure; hold on;

quiver(S, Z, data_ds, data_dz);
title('Direction Field of Humans vs Zombies');
xlabel('Humans');
ylabel('Zombies');
print(fig, '-dpng', 'direction_field');

close(fig);

fig = figure; hold on;
tspan  = [0 35]; % days
y0     = [49999; 1];
[t, y] = ode45(@(t, y) szr(t, y, alpha, beta, gamma, N0), tspan, y0);
plot(t, y(:, 1));
plot(t, y(:, 2));
legend('Human population', 'Zombie population');
title('Human vs Zombie Poulations');
xlabel('Days');
ylabel('Number of individuals');
