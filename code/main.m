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
                    (beta-alpha).*Z - gamma, (beta-alpha).*S - gamma];

%%%
% equilibirum points
% no humans, only zombies
equilibirum_jacobian_1 = Jacobian(0, N0);
e_vals_1 = eig(equilibirum_jacobian_1);

% no zombies, only humans
equilibirum_jacobian_2 = Jacobian(N0, 0);
e_vals_2 = eig(equilibirum_jacobian_2);

%%%
% Stabliity of equilibirum points
disp('Without an antidote:');
if any(e_vals_1 > 0) % all zombies
  disp('The all-zombies equilibirum point is asymptotically unstable');
else
  disp('The all-zombies equilibirum point is asymptotically stable');
end
% Stabliity of equilibirum points
if any(e_vals_2 > 0) % all humans
  disp('The all-humans equilibirum point is asymptotically unstable');
else
  disp('The all-humans equilibirum point is asymptotically stable');
end

%%%
% plot humans and zombies as direction field
[S, Z] = meshgrid(0:(N0/10):N0, 0:(N0/10):N0);
data_ds = dS(S, Z);
data_dz = dZ(S, Z);

fig = figure; hold on;
quiver(S, Z, data_ds, data_dz);
title('Direction Field of Humans vs Zombies');
xlabel('Humans');
ylabel('Zombies');
print(fig, '-dpng', 'direction_field');

%%%
% plot humans and zombies over time
fig = figure; hold on;
tspan  = [0 35]; % days
y0     = [49999; 1];
[t, y] = ode45(@(t, y) szr(t, y, alpha, beta, gamma, N0), tspan, y0);
plot(t, y(:, 1)); % S
plot(t, y(:, 2)); % Z
plot(t, N0 - y(:, 1) - y(:, 2)); % R
legend('Human population', 'Zombie population', 'Removed population');
title('Human vs Zombie Poulations');
xlabel('Days');
ylabel('Number of individuals');
print(fig, '-dpng', 'without_antidote');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SZR Model 2
rho = 0.4;

%%%
% Jacobian with antidote
Jacobian_2 = @(S, Z) [ -beta.*Z,                -beta.*S + rho
                       (beta-alpha).*Z - gamma, (beta-alpha).*S - gamma - rho];

%%%
% equilibirum points
% no zombies, only humans
equilibirum_jacobian_1 = Jacobian_2(N0, 0);
e_vals = eig(equilibirum_jacobian_1);
disp('With an antidote:');
if any(e_vals > 0) % all zombies
  disp('The all-humans equilibirum point is asymptotically unstable');
else
  disp('The all-humans equilibirum point is asymptotically stable');
end
% % alternate point -> NOT IN RANGE - ZOMBIE NUMBER IS NEGATIVE
% equilibirum_jacobian_2 = Jacobian_2(rho/beta, N0*(gamma - rho/beta)/(gamma + (alpha*rho)/beta) );
% e_vals_2 = eig(equilibirum_jacobian_1);


%%%
% plot humans and zombies as direction field
[S, Z]  = meshgrid(0:(N0/10):N0, 0:(N0/10):N0);
dS      = @(S, Z) -beta.*Z.*S + rho.*Z;
dZ      = @(S, Z) (beta - alpha).*Z.*S + gamma.*(N0 - S - Z) - rho.*Z;
data_ds = dS(S, Z);
data_dz = dZ(S, Z);

fig = figure; hold on;
quiver(S, Z, data_ds, data_dz);
title('Direction Field of Humans vs Zombies (with Antidote)');
xlabel('Humans');
ylabel('Zombies');
print(fig, '-dpng', 'direction_field_with_antidote');


%%%
% plot humans and zombies as direction field
fig = figure; hold on;
tspan  = [0 35];                 % days
y0     = [49999; 1];
[t, y] = ode45(@(t, y) szr_with_antidote(t, y, alpha, beta, gamma, rho, N0), tspan, y0);
plot(t, y(:, 1));                % S
plot(t, y(:, 2));                % Z
plot(t, N0 - y(:, 1) - y(:, 2)); % R
legend('Human population', 'Zombie population', 'Removed population');
title('Human vs Zombie Poulations with Antidote');
xlabel('Days');
ylabel('Number of individuals');
print(fig, '-dpng', 'with_antidote');
