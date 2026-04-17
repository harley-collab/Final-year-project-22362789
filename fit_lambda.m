  function lam = fit_lambda(t_canonical, mean_combo, t_start, t_end_win)
    mask = t_canonical >= t_start & t_canonical <= t_end_win & mean_combo > 0 & isfinite(log(mean_combo));
    t_fit = t_canonical(mask);%time window
    y_fit = log(mean_combo(mask));
    t_bar = mean(t_fit);%mean
    y_bar = mean(y_fit);%mean
    S_tt  = sum((t_fit - t_bar).^2);% Variance of time
    S_ty  = sum((t_fit - t_bar) .* (y_fit - y_bar));%covariance 
    lam   = S_ty / S_tt;%LCE value
  end