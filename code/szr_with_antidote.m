function dydt = szr_with_antidote(~, y, alpha, beta, gamma, rho, N0)

  dydt = [-beta * y(1) * y(2) + rho.*y(2);
           beta * y(1) * y(2) + gamma * (N0 - y(1) - y(2)) - alpha * y(1) * y(2) - rho.*y(2);];

end
