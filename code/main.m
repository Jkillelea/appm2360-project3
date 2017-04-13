clear; clc; close all;

% params
alpha = 0.00001;
beta  = 0.00003;
gamma = 0.00001;
N0    = 5000;

% SZR Model
dS = @(S, Z) -beta.*Z.*S;
dZ = @(S, Z) beta.*Z.*S + gamma.*(N0 - S - Z) - alpha.*Z.*S;

data_ds = ones(N0);
data_dz = ones(N0);

for S = 1:length(N0)
  for Z = 1:length(N0)
    data_ds(S, Z) = dS(S, Z);
    data_dz(S, Z) = dZ(S, Z);
  end
end



% tspan  = [0 35]; % days
% y0     = [49999; 1];
% [t, y] = ode45(@(t, y) szr(t, y, alpha, beta, gamma, N0), tspan, y0);
%
% plot(t, y);
