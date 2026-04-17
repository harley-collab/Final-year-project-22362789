function C3_safe = findCriticalC3_NR(theta0, C1, C2, m, R, h, g, I1, I3, phi_dot0, psi_dot0, tspan, r_Gc, y0, tol)

  if nargin < 17
    tol = 1e-6;
  end

  Froot = @(C3) computeFroot(C3, theta0, C1, C2, m, R, h, g, I1, I3, phi_dot0, psi_dot0, tspan, r_Gc, y0);



  C3 = 0.0001; %start
  step = 0.00001;
  Fval = Froot(C3);



  maxC3 = 0.1;
  found = false;

  while C3 < maxC3
    C3_next = C3 + step;
    F_next = Froot(C3_next);

    if isnan(F_next)
      C3_safe = NaN;
      return;
    end

    if Fval > 0 && F_next <= 0
      C3_lo = C3;
      C3_hi = C3_next;
      found = true;
      break;
    end

    C3 = C3_next;
    Fval = F_next;
  end

  if ~found
    C3_safe = NaN;
    return;
  end


  C3 = C3_hi; %top
  %C3 = C3_lo; %bottom
  fd_h = 1e-5;


  

  for iter = 1:100
    F = Froot(C3);
    if isnan(F)
      C3_safe = NaN;
      return;
    end

    if abs(F) < tol
      C3_safe = C3;
      return;
    end

    
   Fp = (Froot(C3 + fd_h) - Froot(C3 - fd_h)) / (2*fd_h);%newton raphson

    if isnan(Fp) || abs(Fp) < 1e-12
      
      C3_safe = C3_hi;%top
      %C3_safe = C3_lo;%bottom
      return;
    end

   
    C3_new = C3 - F/Fp;

   
    if C3_new <= C3_lo || C3_new >= C3_hi
      C3_new = 0.5*(C3 + C3_hi);%top
      %C3_new = 0.5*(C3 + C3_lo);%bottom
    end

    if abs(C3_new - C3) < tol
      C3_safe = C3_new;
      return;
    end

    C3 = C3_new;
  end

  
  C3_safe = C3;
end