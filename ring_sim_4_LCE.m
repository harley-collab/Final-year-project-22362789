function ring_sim_4_LCE()
  tic;
  basePlotDir = 'plots_4_4_1';
  
  if ~exist(basePlotDir,'dir'), mkdir(basePlotDir); end
  

  dataDir = 'data_4_4_1';
  if ~exist(dataDir,'dir'), mkdir(dataDir); end

  pool = gcp('nocreate');%multithread
  if isempty(pool)
    parpool('local', feature('numcores'));
  end

  m  = 1.0;
  R  = 1.025;
  h  = 0.88;
  g  = 9.81;
  I1 = 0.5772;
  I3 = 1.0253;
  r_Gc = [0; -h/2; R];

  Cv = [0.05  0.05	0.01];%drag
  C1 = Cv(1); C2 = Cv(2); C3 = Cv(3);
  tag = sprintf('C1_%g_C2_%g_C3_%g', C1, C2, C3);

  dt    = 0.001;%time step
  t_end = 600;
  nSteps = floor(t_end / dt);%No. steps

  nMC   = 1000;%runcount
  sigma = 1e-10;%perturbation

  theta0_nom   = 0.55;%IC values
  phi_dot0_nom = 10;
  psi_dot0_nom = 0.2;

  rng(42);
  pert3 = sigma * randn(nMC-1, 3);%pertubation

  all_theta0   = [theta0_nom;   theta0_nom   + pert3(:,1)];%IC pert
  all_phi_dot0 = [phi_dot0_nom; phi_dot0_nom + pert3(:,2)];
  all_psi_dot0 = [psi_dot0_nom; psi_dot0_nom + pert3(:,3)];

  all_x      = cell(nMC,1);
  all_y_traj = cell(nMC,1);
  all_z      = cell(nMC,1);
  all_t      = cell(nMC,1);
  all_theta  = cell(nMC,1);
  all_phi    = cell(nMC,1);
  all_psi    = cell(nMC,1);

  dq = parallel.pool.DataQueue;
  completedMC = 0;
  afterEach(dq, @(~) progressCallback());

  function progressCallback()
    completedMC = completedMC + 1;
    pct = completedMC / nMC * 100;
    fprintf('\rProgress: %5.1f%% (%d/%d)', pct, completedMC, nMC);
    if completedMC == nMC
      fprintf('\n');
    end
  end

  parfor iMC = 1:nMC

    theta0   = all_theta0(iMC);
    phi_dot0 = all_phi_dot0(iMC);
    psi_dot0 = all_psi_dot0(iMC);

    w10 = 0;
    w20 = psi_dot0 + phi_dot0 * sin(theta0);
    w30 = phi_dot0 * cos(theta0);

    y = [0; 0; R*cos(theta0); w10; w20; w30; theta0; 0; 0];

    t_vals = zeros(nSteps+1,1);
    y_vals = zeros(nSteps+1,9);
    y_vals(1,:) = y';

    nFinal = nSteps;
%fixed time step variation
    for n = 1:nSteps
        dydt = eom([], y, m, R, h, g, I1, I3, C1, C2, C3, r_Gc);

        v_dot     = dydt(1:3);
        omega_dot = dydt(4:6);

        y(4:6) = y(4:6) + dt * omega_dot;

        w1    = y(4);
        w3    = y(6);
        theta = y(7);

        phi_dot   = w3 / cos(theta);
        theta_dot = w1;
        psi_dot   = y(5) - phi_dot * sin(theta);

        y(7) = y(7) + dt * theta_dot;
        y(8) = y(8) + dt * phi_dot;
        y(9) = y(9) + dt * psi_dot;
        y(1:3) = y(1:3) + dt * v_dot;

        y_vals(n+1,:) = y';
        t_vals(n+1)   = n * dt;

        if y(7) >= pi*9999/20000 || y(7) <= 0 %event
            nFinal = n;
            break;
        end
    end

    all_x{iMC}      = y_vals(1:nFinal+1, 1);
    all_y_traj{iMC} = y_vals(1:nFinal+1, 2);
    all_z{iMC}      = y_vals(1:nFinal+1, 3);
    all_t{iMC}      = t_vals(1:nFinal+1);
    all_theta{iMC}  = y_vals(1:nFinal+1, 7);
    all_phi{iMC}    = y_vals(1:nFinal+1, 8);
    all_psi{iMC}    = y_vals(1:nFinal+1, 9);

    send(dq, iMC);
  end

 
  n_ref = length(all_t{1});%refrence trajectory
  d_combo = NaN(n_ref, nMC-1);

  for k = 2:nMC%compare to ref
    n_k      = length(all_t{k});
    n_common = min(n_ref, n_k);

    dth = all_theta{k}(1:n_common) - all_theta{1}(1:n_common);%diference
    dph = all_phi{k}(1:n_common)   - all_phi{1}(1:n_common);
    dps = all_psi{k}(1:n_common)   - all_psi{1}(1:n_common);

    d_combo(1:n_common, k-1) = sqrt(dth.^2 + dph.^2 + dps.^2);%euclidean distance
  end

  mean_combo = mean(d_combo, 2, 'omitnan');%mean divergence
  t_canonical = all_t{1};

  win_pts   = round(0.02 * length(t_canonical));%for slope
  step_pts  = max(1, round(0.005 * length(t_canonical)));  
  sat_thresh = 0.2;   

  valid_mask = mean_combo > 0 & isfinite(log(mean_combo));
  log_d = NaN(size(mean_combo));
  log_d(valid_mask) = log(mean_combo(valid_mask));

  
  early_slope = NaN;%estimate linear region
  for ii = 1 : step_pts : (length(t_canonical) - win_pts + 1)
    idx = ii : ii+win_pts-1;
    sub_t = t_canonical(idx);
    sub_y = log_d(idx);
    ok = isfinite(sub_y);
    if sum(ok) >= 10
      p_early = polyfit(sub_t(ok), sub_y(ok), 1);
      early_slope = p_early(1);
      break;
    end
  end

  t_sat = NaN;%saturation
  if ~isnan(early_slope) && early_slope > 0
    for ii = 1 : step_pts : (length(t_canonical) - win_pts + 1)
      idx = ii : ii+win_pts-1;
      sub_t = t_canonical(idx);
      sub_y = log_d(idx);
      ok = isfinite(sub_y);
      if sum(ok) < 10, continue; end
      p_local = polyfit(sub_t(ok), sub_y(ok), 1);
      if p_local(1) <= sat_thresh * early_slope
        t_sat = t_canonical(ii + floor(win_pts/2)); 
        break;
      end
    end
  end

  if ~isnan(t_sat)
    fprintf('\n  Saturation time  t_sat = %.1f s\n', t_sat);
  else
    fprintf('\n  Saturation time  : not detected within simulation window\n');
  end

  
  mask   = mean_combo > 0;
  t_fit  = t_canonical(mask);
  y_fit  = log(mean_combo(mask));

  coeffs = polyfit(t_fit, y_fit, 1);
  lambda = coeffs(1);%LCE
  A_fit  = exp(coeffs(2));

  fprintf('\nLambda    = %.5f s^-1\n', lambda);
  fprintf('1/lambda  = %.2f s\n',      1/lambda);

  lambda_plot = fit_lambda(t_canonical, mean_combo, 0, 10);
  mask_plot   = t_canonical >= 0 & t_canonical <= 10 & mean_combo > 0;
  t_plot5     = t_canonical(mask_plot);
  A_plot5     = mean_combo(find(mask_plot, 1)) / exp(lambda_plot * t_plot5(1));

  fig1 = figure('Visible','off','Position',[100 100 1000 560]);
  semilogy(t_canonical, mean_combo, 'b', 'LineWidth', 2); hold on;
  semilogy(t_canonical, A_plot5 * exp(lambda_plot * t_canonical), 'r--', 'LineWidth', 1.5);
  xlabel('Time (s)', 'FontSize', 18);
  ylabel('Mean divergence ||\DeltaIC||', 'FontSize', 18);
  title(sprintf('Lyapunov exponent  \\lambda = %.5f s^{-1} \\sigma = %.0e', lambda, sigma), 'FontSize', 18);
  legend('Mean divergence', ...
         sprintf('Exponential fit %.2e \\cdot exp(t / %.2f)', A_fit, 1/lambda), ...
         'Location','northwest','FontSize',14);
  set(gca, 'FontSize', 16);
  ylim([1e-10, 1e-6]);
  grid on;

  fname1 = fullfile(basePlotDir, ...
    ['traj_2D_XY_Lyapunov_vector_perturbation_tfit_' tag '.png']);

  if exist('exportgraphics','builtin') || exist('exportgraphics','file')
    exportgraphics(fig1, fname1, 'Resolution', 150);
  else
    print(fig1, fname1, '-dpng', '-r150');
  end
  close(fig1);


  
  t_fit_end   = t_end;
  mask   = t_canonical >= t_fit_start & t_canonical <= t_fit_end ...
         & mean_combo > 0 & isfinite(log(mean_combo));

  t_fit = t_canonical(mask);
  y_fit = log(mean_combo(mask));
  n_fit = length(t_fit);

  t_bar = mean(t_fit);
  y_bar = mean(y_fit);
  S_tt  = sum((t_fit - t_bar).^2);
  S_ty  = sum((t_fit - t_bar) .* (y_fit - y_bar));
  S_yy  = sum((y_fit - y_bar).^2);

  lambda_full = S_ty / S_tt;%LCE
  A_fit       = exp(y_bar - lambda_full * t_bar);

  y_hat  = lambda_full * t_fit + log(A_fit);
  SS_res = sum((y_fit - y_hat).^2);
  R2     = 1 - SS_res / S_yy;%eqn. 5.9

  s2           = SS_res / (n_fit - 2);
  sigma_lambda = sqrt(s2 / S_tt);%eqn. 5.10

  t_stat    = lambda_full / sigma_lambda;%eqn. 5.23
  nu        = n_fit - 2;
  t_crit    = tinv(0.95, nu);
  reject_H0 = t_stat > t_crit;%null hypothesis

  t_restrict_end    = t_end / 2;
  lambda_restricted = fit_lambda(t_canonical, mean_combo, 0, t_restrict_end);

  delta_lambda = abs(lambda_full - lambda_restricted);%eqn. 5.25

  lambda = lambda_full;

  fprintf('    n_fit        : %d points                 \n', n_fit);
  fprintf('    lambda       : %+.5f s^-1               \n', lambda);
  fprintf('    lambda_full  : %+.5f s^-1               \n', lambda_full);
  fprintf('    lambda_restr : %+.5f s^-1               \n', lambda_restricted);
  fprintf('    Delta lambda : %.5f s^-1               \n', delta_lambda);
  fprintf('    sigma_lambda : %.5f s^-1                \n', sigma_lambda);
  fprintf('    R^2          : %.6f                     \n', R2);
  fprintf('    t-statistic  : %.3f                     \n', t_stat);
  fprintf('    t_crit(0.95) : %.3f                     \n', t_crit);
  fprintf('    Reject H0    : %s                        \n', string(reject_H0));
  fprintf('    A (prefactor): %.3e                     \n', A_fit);
  fprintf('    1/lambda     : %.2f s                   \n', 1/lambda);
  if ~isnan(t_sat)
    fprintf('    t_sat        : %.1f s                  \n', t_sat);
  else
    fprintf('    t_sat        : not detected             \n');
  end
 

  fprintf('Plot saved: %s\n', fname1);
  fprintf('\nDONE (Euler-Cromer)  elapsed = %.1f s\n', toc);

end
