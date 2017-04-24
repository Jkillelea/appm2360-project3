clear; clc; close all;

% params
alpha = 0.00001;
beta  = 0.00003;
gamma = 0.00001;
N0    = 50000;

% SZR Model 1
dS = @(S, Z) -beta.*Z.*S;
dZ = @(S, Z) (beta - alpha).*Z.*S + gamma.*(N0 - S - Z);

% jacobian, evaluated at both equilibirum points
Jacobian = @(S, Z) [-beta.*Z,                          -beta.*S;
                    (beta-alpha).*Z + gamma.*(N0 - Z), (beta-alpha).*S + gamma.*(N0 - S)];

% no humans, only zombies
equilibirum_jacobian_1 = Jacobian(0, N0);
% no zombies, only humans
equilibirum_jacobian_2 = Jacobian(N0, 0);


[S, Z] = meshgrid(0:(N0/10):N0, 0:(N0/10):N0);

data_ds = dS(S, Z);
data_dz = dZ(S, Z);

fig = figure; hold on;
quiver(S, Z, data_ds, data_dz);
title('Direction Field of Humans vs Zombies');
xlabel('Humans');
ylabel('Zombies');
print(fig, '-dpng', 'direction_field');

% close(fig);

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
print(fig, '-dpng', 'without_antidote');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SZR Model 2
rho = 0.4;
dS = @(S, Z) -beta.*Z.*S + rho.*Z;
dZ = @(S, Z) (beta - alpha).*Z.*S + gamma.*(N0 - S - Z) - rho.*Z;
% plug in and solve with ode45 again

fig = figure; hold on;
tspan  = [0 35]; % days
y0     = [49999; 1];
[t, y] = ode45(@(t, y) szr_with_antidote(t, y, alpha, beta, gamma, rho, N0), tspan, y0);
plot(t, y(:, 1));
plot(t, y(:, 2));
legend('Human population', 'Zombie population');
title('Human vs Zombie Poulations with Antidote');
xlabel('Days');
ylabel('Number of individuals');
print(fig, '-dpng', 'with_antidote');
