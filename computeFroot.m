function Fval = computeFroot(C3, theta0, C1, C2, m, R, h, g, I1, I3, phi_dot0, psi_dot0, tspan, r_Gc, y0)
%back up root finder
  try
    [~, ~, Ff_traj] = runSimulation(theta0, C1, C2, C3, m, R, h, g, I1, I3, phi_dot0, psi_dot0, tspan, r_Gc, y0);
    Fval = min(Ff_traj);
  catch
    Fval = NaN;
  end
end