function ring_sim_4_TOT()

tic;
  basePlotDir = 'plots_4_2_1';
  dataDir     = 'data_4_2_1';

  if ~exist(basePlotDir,'dir'), mkdir(basePlotDir); end
  if ~exist(dataDir,'dir'),     mkdir(dataDir);     end

C1_values = [0.01 0.025 0.05 0.075 0.1 0.25 0.5];%drag range
C2_values = [0.01 0.025 0.05 0.075 0.1 0.25 0.5];
C3_values = [0.01 0.025 0.05 0.075 0.1 0.25 0.5];



m = 1.0;
R = 1.025;
h = 0.88;
g = 9.81;
I1 = 0.5772;
I3 = 1.0253;
r_Gc = [0; -h/2; R];
theta0     = 0.55;
phi0       = 0;
psi0       = 0;
theta_dot0 = 0;
phi_dot0   = 20;
psi_dot0   = 0.5;

tspan = [0 600];

options = odeset('RelTol',1e-5,'AbsTol',1e-7,'MaxStep',1e-1);

C1_all = [];
C2_all = [];
C3_all = [];
retro_all = [];
t_trans_all = [];

totalRuns = numel(C1_values)*numel(C2_values)*numel(C3_values);
runCount = 0;

for C1 = C1_values
for C2 = C2_values
for C3 = C3_values

runCount = runCount + 1;
fprintf('Run %d / %d (%.1f%%)\n', runCount, totalRuns,100*runCount/totalRuns);

w10 = 0;
w20 = psi_dot0 + phi_dot0*sin(theta0);
w30 = phi_dot0*cos(theta0);

y0 = [0;0;R*cos(theta0); w10;w20;w30; theta0;phi0;psi0];

[t,y] = ode45(@(t,y) eom(t,y, m, R, h, g, I1, I3, C1, C2, C3, r_Gc),tspan,y0,options);

Ff  = zeros(length(t),1);
rcdd = zeros(length(t),3);

for k = 1:length(t)
[Ff(k), rcdd(k,:)] = turntrack(y(k,:),m,R,h,g,I1,I3,C1,C2,C3,r_Gc);
end


idx = find(Ff(1:end-1).*Ff(2:end) < 0,1);

if isempty(idx)

transition_found = false;
t_transition = NaN;

fprintf('%.3f %.3f %.3f : No transition\n',C1,C2,C3);

else

transition_found = true;

t_transition = t(idx) - Ff(idx)*(t(idx+1)-t(idx))/(Ff(idx+1)-Ff(idx));%eqn. 5.3

fprintf('%.3f %.3f %.3f %.6f\n',C1,C2,C3,t_transition);

end


C1_all(end+1) = C1;
C2_all(end+1) = C2;
C3_all(end+1) = C3;
retro_all(end+1) = transition_found;
t_trans_all(end+1) = t_transition;



end
end
end
save(fullfile(dataDir, 'sweep_results.mat'), ...
     'C1_all','C2_all','C3_all','retro_all','t_trans_all');

tlog_all = log10(t_trans_all);

fig = figure('visible','off');
scatter3(C1_all, C2_all, C3_all, 80, tlog_all, "filled")
xlabel("C1"); ylabel("C2"); zlabel("C3")
title("Transition time across parameter space")
colormap(jet); colorbar; grid on; view(40,25)
print(fig, fullfile(basePlotDir, 'tot_surf.png'), '-dpng', '-r150');
close(fig);


C3_vals = unique(C3_all);
max_t = zeros(length(C3_vals),1);
min_t = zeros(length(C3_vals),1);
max_C1 = zeros(length(C3_vals),1); max_C2 = zeros(length(C3_vals),1);
min_C1 = zeros(length(C3_vals),1); min_C2 = zeros(length(C3_vals),1);

for i = 1:length(C3_vals)
    c3 = C3_vals(i);
    idx = (C3_all == c3);

    C1s   = C1_all(idx);
    C2s   = C2_all(idx);
    ts    = tlog_all(idx);
    tvals = t_trans_all(idx);

    if sum(idx) < 3; continue; end

    [max_t(i), max_idx] = max(tvals);
    max_C1(i) = C1s(max_idx); max_C2(i) = C2s(max_idx);
    [min_t(i), min_idx] = min(tvals);
    min_C1(i) = C1s(min_idx); min_C2(i) = C2s(min_idx);

    n = 60;
    xlin = linspace(min(C1s), max(C1s), n);
    ylin = linspace(min(C2s), max(C2s), n);
    [X,Y] = meshgrid(xlin, ylin);
    Z = griddata(C1s, C2s, ts, X, Y, "linear");

    fig = figure('visible','off');
    surf(X,Y,Z)
    xlabel("C1"); ylabel("C2"); zlabel("log10(t)")
    title(sprintf("Surface (C3 = %.3f)", c3))
    shading interp; colormap(jet); colorbar; grid on
    filename = sprintf('tot_surf_lin_C3_%0.3f.png', c3);
    print(fig, fullfile(basePlotDir, filename), '-dpng', '-r150');

    C1log = log10(C1s); C2log = log10(C2s);
    xlin = linspace(min(C1log), max(C1log), n);
    ylin = linspace(min(C2log), max(C2log), n);
    [X,Y] = meshgrid(xlin, ylin);
    Z = griddata(C1log, C2log, ts, X, Y, "linear");

    fig = figure('visible','off');
    surf(X,Y,Z)
    xlabel("log10(C1)"); ylabel("log10(C2)"); zlabel("log10(t)")
    title(sprintf("Log Surface (C3 = %.3f)", c3))
    shading interp; colormap(jet); colorbar; grid on
    filename = sprintf('tot_surf_log_C3_%0.3f.png', c3);
    print(fig, fullfile(basePlotDir, filename), '-dpng', '-r150');
end

fig = figure('visible','off');
plot(C3_vals, max_t, '-or', 'LineWidth', 2); hold on
plot(C3_vals, min_t, '-ob', 'LineWidth', 2)
xlabel("C3"); ylabel("Transition time")
title("Max and Min Transition Time vs C3")
legend("Max","Min"); grid on
print(fig, fullfile(basePlotDir, 'tot_max_min.png'), '-dpng', '-r150');
close(fig);
fprintf('\nCompleted %d runs in %.2f seconds\n',runCount,toc);

end