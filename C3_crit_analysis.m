function C3_crit_analysis() 

bottom_files = { 'C3_crit_matrix_bottom.dat', ... % 'bottom' data files
                 'C3_crit_matrix_bottom_2.dat', ...
                 'C3_crit_matrix_bottom_3.dat', ...
                 'C3_crit_matrix_bottom_4.dat', ...
                 'C3_crit_matrix_bottom_5.dat', ...
                 'C3_crit_matrix_bottom_6.dat' };

top_files    = { 'C3_crit_matrix_top.dat', ... % 'top' 
                 'C3_crit_matrix_top_2.dat', ...
                 'C3_crit_matrix_top_3.dat', ...
                 'C3_crit_matrix_top_4.dat', ...
                 'C3_crit_matrix_top_5.dat', ...
                 'C3_crit_matrix_top_6.dat' };

data_bottom = load_all(bottom_files); % Load and merge
data_top    = load_all(top_files);    
    
delta = data_top(:,3) - data_bottom(:,3); % eqn. 5.20

mean_delta  = mean(delta); %eqn. 5.21
sigma_delta = std(delta);  %eqn. 5.22

fprintf('mean(dC3)  = %.6e\n', mean_delta); 
fprintf('sigma(dC3) = %.6e\n', sigma_delta); 
fprintf('Loaded %d bottom points, %d top points\n', size(data_bottom,1), size(data_top,1));

[min_delta, min_idx] = min(abs(delta)); 

C1_min = data_bottom(min_idx, 1); 
C2_min = data_bottom(min_idx, 2); 
C3_bot = data_bottom(min_idx, 3); 
C3_top = data_top(min_idx, 3);    

fprintf('Min |dC3|        = %.16e\n', min_delta); 
fprintf('Location (C1,C2) = (%.6f, %.6f)\n', C1_min, C2_min); 
fprintf('C3_crit bottom   = %.6e\n', C3_bot); 
fprintf('C3_crit top      = %.6e\n', C3_top); 

fprintf('\n BOTTOM surface\n'); 
[beta_b, se_beta_b, gam_b, se_gam_b, A_b, R2_b] = fit_powerlaw(data_bottom); % Fit data

fprintf('\n TOP surface \n'); 
[beta_t, se_beta_t, gam_t, se_gam_t, A_t, R2_t] = fit_powerlaw(data_top); 


fprintf('%-10s  %-20s  %-20s  %-10s  %-10s\n', ... 
        'Dataset', 'beta', 'gamma', 'A', 'R^2');
fprintf('%-10s  %-20s  %-20s  %-10s  %-10s\n', ... 
        '-------', '----', '-----', '-', '---');

fprintf('%-10s  %+.4f +/- %.4f   %+.4f +/- %.4f   %.4f   %.6f\n','Bottom', beta_b, se_beta_b, gam_b, se_gam_b, A_b, R2_b);

fprintf('%-10s  %+.4f +/- %.4f   %+.4f +/- %.4f   %.4f   %.6f\n','Top',    beta_t, se_beta_t, gam_t, se_gam_t, A_t, R2_t);

end 


function data = load_all(filelist) % Function to load multiple files
    data = []; 
    for k = 1:length(filelist) 
        raw = importdata(filelist{k}, ' ', 2); % header lines
        data = [data; raw.data]; 
    end
end 


function [beta, se_beta, gam, se_gam, A, R2] = fit_powerlaw(data) 

    C1    = data(:,1); 
    C2    = data(:,2); 
    C3c   = data(:,3); 
    n     = length(C3c); 

    y  = log(C3c); % Log for linear version of fit for eqn. 5.16
    X  = [ones(n,1), log(C1), log(C2)]; 

    c    = X \ y; 
    lnA  = c(1);  %eqn. 5.19
    beta = c(2);  % eqn. 5.17
    gam  = c(3);  % eqn. 5.18
    A    = exp(lnA); % eqn. 5.19

    yhat   = X * c; 
    resid  = y - yhat; % Residuals
    SS_res = sum(resid.^2); % Residual sum of squares
    SS_tot = sum((y - mean(y)).^2); % Total sum of squares
    R2     = 1 - SS_res / SS_tot; %eqn. 5.9

    s2      = SS_res / (n - 3); % Estimate of variance (n - parameters)
    cov_c   = s2 * inv(X'*X); % Covariance matrix of coefficients
    se_beta = sqrt(cov_c(2,2)); % eqn. 5.10
    se_gam  = sqrt(cov_c(3,3)); 

    fprintf('  beta  = %+.6f +/- %.6f\n', beta,  se_beta); 
    fprintf('  gamma = %+.6f +/- %.6f\n', gam,   se_gam); 
    fprintf('  A     = %.6f\n',            A); 
    fprintf('  R^2   = %.6f\n',            R2); 
    fprintf('  n     = %d\n',              n); 

end 