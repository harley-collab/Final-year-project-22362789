function [t,y,Ff] = runSimulation(theta0, C1, C2, C3, m, R, h, g, I1, I3, phi_dot0, psi_dot0, tspan, r_Gc, y0)

  teval = (tspan(1) : 0.5 : tspan(2))';

  options = odeset('RelTol',1e-5,'AbsTol',1e-5,'MaxStep',1e-1, 'Events', @(t,y) combinedEvents(t,y) );

  [t,y] = ode45(@(t,y) eom(t,y,m,R,h,g,I1,I3,C1,C2,C3,r_Gc), teval, y0, options);

  Ff = zeros(size(t));
  for i = 1:length(t)
    Ff(i) = turntrack(y(i,:),m,R,h,g,I1,I3,C1,C2,C3,r_Gc);
  end
end