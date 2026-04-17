% for producing plots of phase change boundary and picking points

files_top = {
    'C:/Users/james/Documents/MATLAB/C3_crit_matrix_top.dat'
    'C:/Users/james/Documents/MATLAB/C3_crit_matrix_top_2.dat'
    'C:/Users/james/Documents/MATLAB/C3_crit_matrix_top_3.dat'
    'C:/Users/james/Documents/MATLAB/C3_crit_matrix_top_4.dat'
    'C:/Users/james/Documents/MATLAB/C3_crit_matrix_top_5.dat'
    'C:/Users/james/Documents/MATLAB/C3_crit_matrix_top_6.dat'
};

files_bottom = {
    'C:/Users/james/Documents/MATLAB/C3_crit_matrix_bottom.dat'
    'C:/Users/james/Documents/MATLAB/C3_crit_matrix_bottom_2.dat'
    'C:/Users/james/Documents/MATLAB/C3_crit_matrix_bottom_3.dat'
    'C:/Users/james/Documents/MATLAB/C3_crit_matrix_bottom_4.dat'
    'C:/Users/james/Documents/MATLAB/C3_crit_matrix_bottom_5.dat'
    'C:/Users/james/Documents/MATLAB/C3_crit_matrix_bottom_6.dat'
};

data_top    = loadAllFiles(files_top);
data_bottom = loadAllFiles(files_bottom);

C1_values = unique(data_top(:,1));
C2_values = unique(data_top(:,2));
[C2_grid, C1_grid] = meshgrid(C2_values, C1_values);

C3_top    = loadMatrix_fromData(data_top);
C3_bottom = loadMatrix_fromData(data_bottom);
C3_diff   = C3_top - C3_bottom;

C3_pos = C3_diff; C3_pos(C3_diff < 0) = 0;
C3_neg = C3_diff; C3_neg(C3_diff > 0) = 0;
Z_zero = zeros(size(C3_diff));

nC1 = size(C3_diff, 1);
nC2 = size(C3_diff, 2);
edges = { 1,    1:nC2; nC1,  1:nC2; 1:nC1, 1;1:nC1, nC2};

nLevels = 100;   


query = [0.750, 0.7500];   


[~, iC1] = min(abs(C1_values - query(1)));
[~, iC2] = min(abs(C2_values - query(2)));

C3_top_val    = C3_top(iC1, iC2);
C3_bottom_val = C3_bottom(iC1, iC2);
C3_diff_val   = C3_top_val - C3_bottom_val;

fprintf('\n── Table lookup: C1 = %.6f, C2 = %.6f ──\n', C1_values(iC1), C2_values(iC2));
fprintf('  C3_top    (prograde)   = %.6f\n', C3_top_val);
fprintf('  C3_bottom (retrograde) = %.6f\n', C3_bottom_val);
fprintf('  DeltaC3                = %.6f\n', C3_diff_val);

figure(1)
surf(C2_grid, C1_grid, C3_top)
xlabel('C2'); ylabel('C1'); zlabel('C3_{crit}')
title('Top Surface')
colormap(gca, jet(256)); colorbar; shading interp; grid on; view(45,30)
saveas(figure(1), 'fig1_top_surface.png')

figure(2)
surf(C2_grid, C1_grid, C3_bottom)
xlabel('C2'); ylabel('C1'); zlabel('C3_{crit}')
title('Bottom Surface')
colormap(gca, jet(256)); colorbar; shading interp; grid on; view(45,30)
saveas(figure(2), 'fig2_bottom_surface.png')

figure(3)
h_top = surf(C2_grid, C1_grid, C3_top);
hold on
h_bot = surf(C2_grid, C1_grid, C3_bottom);
hold off
h_top.FaceColor = [0.2 0.6 1.0]; h_top.FaceAlpha = 0.7; h_top.EdgeColor = 'none';
h_bot.FaceColor = [1.0 0.4 0.2]; h_bot.FaceAlpha = 0.7; h_bot.EdgeColor = 'none';
xlabel('C2'); ylabel('C1'); zlabel('C3_{crit}')
title('Top & Bottom Overlay')
legend([h_top h_bot], 'Top','Bottom')
grid on; view(45,30)
saveas(figure(3), 'fig3_overlay.png')

figure(4)
surf(C2_grid, C1_grid, C3_diff)
xlabel('C2'); ylabel('C1'); zlabel('\DeltaC3_{crit}  (Top - Bottom)')
title('Difference: Top minus Bottom')
colormap(gca, jet(256)); colorbar; shading interp; grid on; view(45,30)
hold on
contour3(C2_grid, C1_grid, C3_diff, [0 0], 'k-', 'LineWidth', 2)
hold off
legend('C3 difference', 'Zero crossing', 'Location', 'best')
saveas(figure(4), 'fig4_difference.png')

figure(5); clf
surf(C2_grid, C1_grid, C3_pos, 'FaceColor',[0.2 0.6 1.0],'FaceAlpha',0.6,'EdgeColor','none')
hold on
surf(C2_grid, C1_grid, Z_zero, 'FaceColor',[0.2 0.6 1.0],'FaceAlpha',0.2,'EdgeColor','none')
surf(C2_grid, C1_grid, C3_neg, 'FaceColor',[1.0 0.4 0.2],'FaceAlpha',0.6,'EdgeColor','none')
surf(C2_grid, C1_grid, Z_zero, 'FaceColor',[1.0 0.4 0.2],'FaceAlpha',0.2,'EdgeColor','none')
contour3(C2_grid, C1_grid, C3_diff, [0 0], 'k-', 'LineWidth', 2.5)
for e = 1:4
    iC1 = edges{e,1}; iC2 = edges{e,2};
    xv = C2_grid(iC1,iC2); xv=xv(:)';
    yv = C1_grid(iC1,iC2); yv=yv(:)';
    zv = C3_diff(iC1,iC2); zv=zv(:)';
    zv_pos=zv; zv_pos(zv<0)=0; drawCurtain(xv,yv,zv_pos,[0.2 0.6 1.0],0.3)
    zv_neg=zv; zv_neg(zv>0)=0; drawCurtain(xv,yv,zv_neg,[1.0 0.4 0.2],0.3)
end
hold off
xlabel('C2'); ylabel('C1'); zlabel('\DeltaC3_{crit}  (Top - Bottom)')
title('Difference: Top minus Bottom (shaded volume)')
grid on; view(0,0)
h1=line(NaN,NaN,'Color',[0.2 0.6 1.0],'LineWidth',6,'DisplayName','Top > Bottom');
h2=line(NaN,NaN,'Color',[1.0 0.4 0.2],'LineWidth',6,'DisplayName','Bottom > Top');
h3=line(NaN,NaN,'Color','k','LineWidth',2.5,'DisplayName','Zero crossing');
legend([h1 h2 h3],'Location','best')
saveas(figure(5), 'fig5_shaded_volume.png')

figure(6)
surf(C2_grid, C1_grid, log10(C3_top))
set(gca, 'XScale','log', 'YScale','log')
xlabel('C2'); ylabel('C1'); zlabel('log_{10}(C3_{crit})')
title('Top Surface (log axes)')
colormap(gca, jet(256)); colorbar; shading interp; grid on; view(45,30)
saveas(figure(6), 'fig6_top_log.png')

figure(7)
surf(C2_grid, C1_grid, log10(C3_bottom))
set(gca, 'XScale','log', 'YScale','log')
xlabel('C2'); ylabel('C1'); zlabel('log_{10}(C3_{crit})')
title('Bottom Surface (log axes)')
colormap(gca, jet(256)); colorbar; shading interp; grid on; view(45,30)
saveas(figure(7), 'fig7_bottom_log.png')

figure(8)
h_top2 = surf(C2_grid, C1_grid, log10(C3_top));
hold on
h_bot2 = surf(C2_grid, C1_grid, log10(C3_bottom));
hold off
set(gca, 'XScale','log', 'YScale','log')
h_top2.FaceColor = [0.2 0.6 1.0]; h_top2.FaceAlpha = 0.7; h_top2.EdgeColor = 'none';
h_bot2.FaceColor = [1.0 0.4 0.2]; h_bot2.FaceAlpha = 0.7; h_bot2.EdgeColor = 'none';
xlabel('C2'); ylabel('C1'); zlabel('log_{10}(C3_{crit})')
title('Top & Bottom Overlay (log axes)')
legend([h_top2 h_bot2], 'Top','Bottom')
grid on; view(45,30)
saveas(figure(8), 'fig8_overlay_log.png')

figure(9)
surf(C2_grid, C1_grid, C3_diff)
set(gca, 'XScale','log', 'YScale','log')
xlabel('C2'); ylabel('C1'); zlabel('\DeltaC3_{crit}  (Top - Bottom)')
title('Difference: Top minus Bottom (log axes)')
colormap(gca, jet(256)); colorbar; shading interp; grid on; view(45,30)
hold on
contour3(C2_grid, C1_grid, C3_diff, [0 0], 'k-', 'LineWidth', 2)
hold off
legend('C3 difference', 'Zero crossing', 'Location', 'best')
saveas(figure(9), 'fig9_difference_log.png')

figure(10); clf
surf(C2_grid, C1_grid, C3_pos, 'FaceColor',[0.2 0.6 1.0],'FaceAlpha',0.6,'EdgeColor','none')
hold on
surf(C2_grid, C1_grid, Z_zero, 'FaceColor',[0.2 0.6 1.0],'FaceAlpha',0.2,'EdgeColor','none')
surf(C2_grid, C1_grid, C3_neg, 'FaceColor',[1.0 0.4 0.2],'FaceAlpha',0.6,'EdgeColor','none')
surf(C2_grid, C1_grid, Z_zero, 'FaceColor',[1.0 0.4 0.2],'FaceAlpha',0.2,'EdgeColor','none')
contour3(C2_grid, C1_grid, C3_diff, [0 0], 'k-', 'LineWidth', 2.5)
for e = 1:4
    iC1 = edges{e,1}; iC2 = edges{e,2};
    xv = C2_grid(iC1,iC2); xv=xv(:)';
    yv = C1_grid(iC1,iC2); yv=yv(:)';
    zv = C3_diff(iC1,iC2); zv=zv(:)';
    zv_pos=zv; zv_pos(zv<0)=0; drawCurtain(xv,yv,zv_pos,[0.2 0.6 1.0],0.3)
    zv_neg=zv; zv_neg(zv>0)=0; drawCurtain(xv,yv,zv_neg,[1.0 0.4 0.2],0.3)
end
hold off
set(gca, 'XScale','log', 'YScale','log')
xlabel('C2'); ylabel('C1'); zlabel('\DeltaC3_{crit}')
title('Difference: Shaded Volume (log axes)')
grid on; view(45,20)
h1=line(NaN,NaN,'Color',[0.2 0.6 1.0],'LineWidth',6,'DisplayName','Top > Bottom');
h2=line(NaN,NaN,'Color',[1.0 0.4 0.2],'LineWidth',6,'DisplayName','Bottom > Top');
h3=line(NaN,NaN,'Color','k','LineWidth',2.5,'DisplayName','Zero crossing');
legend([h1 h2 h3],'Location','best')
saveas(figure(10), 'fig10_shaded_volume_log.png')

figure(11)
contourf(C1_grid, C2_grid, C3_top, nLevels)
hold on
contour(C1_grid, C2_grid, C3_top, nLevels, 'k-', 'LineWidth', 0.3)
hold off
xlabel('C1'); ylabel('C2')
title('Top Surface (contour)')
colormap(gca, jet(256)); colorbar; axis tight
saveas(figure(11), 'fig11_top_contour.png')


figure(12)
contourf(C1_grid, C2_grid, C3_bottom, nLevels)
hold on
contour(C1_grid, C2_grid, C3_bottom, nLevels, 'k-', 'LineWidth', 0.3)
hold off
xlabel('C1'); ylabel('C2')
title('Bottom Surface (contour)')
colormap(gca, jet(256)); colorbar; axis tight
saveas(figure(12), 'fig12_bottom_contour.png')


figure(13)
contourf(C1_grid, C2_grid, C3_top, nLevels, 'FaceAlpha', 0.6)
hold on
contour(C1_grid, C2_grid, C3_top,    nLevels, 'b-', 'LineWidth', 0.5)
contour(C1_grid, C2_grid, C3_bottom, nLevels, 'r-', 'LineWidth', 0.5)
hold off
xlabel('C1'); ylabel('C2')
title('Top & Bottom Overlay (contour)')
colormap(gca, jet(256)); colorbar; axis tight
h1=line(NaN,NaN,'Color','b','LineWidth',1.5,'DisplayName','Top');
h2=line(NaN,NaN,'Color','r','LineWidth',1.5,'DisplayName','Bottom');
legend([h1 h2],'Location','best')
saveas(figure(13), 'fig13_overlay_contour.png')

figure(14)
contourf(C1_grid, C2_grid, C3_diff, nLevels)
hold on
contour(C1_grid, C2_grid, C3_diff, nLevels, 'k-', 'LineWidth', 0.3)
contour(C1_grid, C2_grid, C3_diff, [0 0],   'w-', 'LineWidth', 2)
hold off
xlabel('C1'); ylabel('C2')
title('Difference: Top minus Bottom (contour)')
colormap(gca, jet(256)); colorbar; axis tight
h1=line(NaN,NaN,'Color','w','LineWidth',2,'DisplayName','Zero crossing');
legend(h1,'Location','best')
saveas(figure(14), 'fig14_difference_contour.png')

figure(15); clf
contourf(C1_grid, C2_grid, C3_pos, nLevels, 'FaceAlpha', 0.8)
colormap(gca, jet(256));
hold on
contourf(C1_grid, C2_grid, C3_neg, nLevels, 'FaceAlpha', 0.8)
contour(C1_grid,  C2_grid, C3_diff, [0 0],  'w-', 'LineWidth', 2.5)
hold off
xlabel('C1'); ylabel('C2')
title('Difference: Shaded Volume (contour)')
colorbar; axis tight
h1=line(NaN,NaN,'Color','w','LineWidth',2.5,'DisplayName','Zero crossing');
legend(h1,'Location','best')
saveas(figure(15), 'fig15_difference_shaded_contour.png')

figure(16)
contourf(C1_grid, C2_grid, log10(C3_top), nLevels)
hold on
contour(C1_grid, C2_grid, log10(C3_top), nLevels, 'k-', 'LineWidth', 0.3)
hold off
set(gca, 'XScale','log', 'YScale','log')
xlabel('C1'); ylabel('C2')
title('Top Surface (log contour)')
colormap(gca, jet(256)); colorbar; axis tight
saveas(figure(16), 'fig16_top_log_contour.png')

figure(17)
contourf(C1_grid, C2_grid, log10(C3_bottom), nLevels)
hold on
contour(C1_grid, C2_grid, log10(C3_bottom), nLevels, 'k-', 'LineWidth', 0.3)
hold off
set(gca, 'XScale','log', 'YScale','log')
xlabel('C1'); ylabel('C2')
title('Bottom Surface (log contour)')
colormap(gca, jet(256)); colorbar; axis tight
saveas(figure(17), 'fig17_bottom_log_contour.png')

figure(18)
contourf(C1_grid, C2_grid, log10(C3_top), nLevels, 'FaceAlpha', 0.6)
hold on
contour(C1_grid, C2_grid, log10(C3_top),    nLevels, 'b-', 'LineWidth', 0.5)
contour(C1_grid, C2_grid, log10(C3_bottom), nLevels, 'r-', 'LineWidth', 0.5)
hold off
set(gca, 'XScale','log', 'YScale','log')
xlabel('C1'); ylabel('C2')
title('Top & Bottom Overlay (log contour)')
colormap(gca, jet(256)); colorbar; axis tight
h1=line(NaN,NaN,'Color','b','LineWidth',1.5,'DisplayName','Top');
h2=line(NaN,NaN,'Color','r','LineWidth',1.5,'DisplayName','Bottom');
legend([h1 h2],'Location','best')
saveas(figure(18), 'fig18_overlay_log_contour.png')

figure(19)
contourf(C1_grid, C2_grid, C3_diff, nLevels)
hold on
contour(C1_grid, C2_grid, C3_diff, nLevels, 'k-', 'LineWidth', 0.3)
contour(C1_grid, C2_grid, C3_diff, [0 0],   'w-', 'LineWidth', 2)
hold off
set(gca, 'XScale','log', 'YScale','log')
xlabel('C1'); ylabel('C2')
title('Difference: Top minus Bottom (log contour)')
colormap(gca, jet(256)); colorbar; axis tight
h1=line(NaN,NaN,'Color','w','LineWidth',2,'DisplayName','Zero crossing');
legend(h1,'Location','best')
saveas(figure(19), 'fig19_difference_log_contour.png')

figure(20); clf
contourf(C1_grid, C2_grid, C3_pos, nLevels, 'FaceAlpha', 0.8)
colormap(gca, jet(256));
hold on
contourf(C1_grid, C2_grid, C3_neg, nLevels, 'FaceAlpha', 0.8)
contour(C1_grid,  C2_grid, C3_diff, [0 0],  'w-', 'LineWidth', 2.5)
hold off
set(gca, 'XScale','log', 'YScale','log')
xlabel('C1'); ylabel('C2')
title('Difference: Shaded Volume (log contour)')
colorbar; axis tight
h1=line(NaN,NaN,'Color','w','LineWidth',2.5,'DisplayName','Zero crossing');
legend(h1,'Location','best')
saveas(figure(20), 'fig20_difference_shaded_log_contour.png')
%picking points in transition region for LCE
C1_vec   = C1_grid(:);
C2_vec   = C2_grid(:);
diff_vec = C3_diff(:);
top_vec  = C3_top(:);
bot_vec  = C3_bottom(:);

nPts = 5;
tol_zero    = 0.05 * (max(diff_vec) - min(diff_vec));  
tol_pos     = 0.05 * max(diff_vec);                    
tol_neg     = 0.05 * abs(min(diff_vec));               
half_pos    = max(diff_vec) / 2;
half_neg    = min(diff_vec) / 2;

function idx = spreadPoints(C1_vec, C2_vec, diff_vec, target, nPts, tol)
    
    candidates = find(abs(diff_vec - target) < tol);
    
    if length(candidates) < nPts
        [~, sorted] = sort(abs(diff_vec - target));
        candidates = sorted(1:min(10*nPts, length(sorted)));
    end

    C1_c = C1_vec(candidates);
    C2_c = C2_vec(candidates);
    
    C1_edges = linspace(min(C1_c), max(C1_c) + eps, nPts+1);
    C2_edges = linspace(min(C2_c), max(C2_c) + eps, nPts+1);
    
    idx = zeros(nPts, 1);
    used = false(length(candidates), 1);
    
    for b = 1:nPts
        in_bin = ~used & ...
                 C1_c >= C1_edges(b) & C1_c < C1_edges(b+1) & ...
                 C2_c >= C2_edges(b) & C2_c < C2_edges(b+1);
        if any(in_bin)
            bin_idx = find(in_bin);
            [~, best] = min(abs(diff_vec(candidates(bin_idx)) - target));
            chosen = bin_idx(best);
        else
            tmp = abs(diff_vec(candidates) - target);
            tmp(used) = inf;
            [~, chosen] = min(tmp);
        end
        idx(b) = candidates(chosen);
        used(chosen) = true;
    end
end

idx_zero    = spreadPoints(C1_vec, C2_vec, diff_vec, 0,        nPts, tol_zero);
idx_pos_mid = spreadPoints(C1_vec, C2_vec, diff_vec, half_pos, nPts, tol_pos);
idx_neg_mid = spreadPoints(C1_vec, C2_vec, diff_vec, half_neg, nPts, tol_neg);



fprintf('\n[1] Zero crossing (C3_top = C3_bottom) — 5 spread points:\n')
fprintf('    %-12s %-12s %-12s %-12s %-12s\n', 'C1','C2','C3_top','C3_bot','C3_diff')
for i = 1:nPts
    fprintf('    %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f\n', ...
        C1_vec(idx_zero(i)), C2_vec(idx_zero(i)), ...
        top_vec(idx_zero(i)), bot_vec(idx_zero(i)), diff_vec(idx_zero(i)))
end

fprintf('\n[2] Halfway between zero and max diff (Top > Bottom) — 5 spread points:\n')
fprintf('    %-12s %-12s %-12s %-12s %-12s\n', 'C1','C2','C3_top','C3_bot','C3_diff')
for i = 1:nPts
    fprintf('    %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f\n', ...
        C1_vec(idx_pos_mid(i)), C2_vec(idx_pos_mid(i)), ...
        top_vec(idx_pos_mid(i)), bot_vec(idx_pos_mid(i)), diff_vec(idx_pos_mid(i)))
end

fprintf('\n[3] Halfway between zero and min diff (Bottom > Top) — 5 spread points:\n')
fprintf('    %-12s %-12s %-12s %-12s %-12s\n', 'C1','C2','C3_top','C3_bot','C3_diff')
for i = 1:nPts
    fprintf('    %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f\n', ...
        C1_vec(idx_neg_mid(i)), C2_vec(idx_neg_mid(i)), ...
        top_vec(idx_neg_mid(i)), bot_vec(idx_neg_mid(i)), diff_vec(idx_neg_mid(i)))
end

function drawCurtain(xv, yv, zv, col, alph)
    for i = 1:length(xv)-1
        xp = [xv(i) xv(i+1) xv(i+1) xv(i)];
        yp = [yv(i) yv(i+1) yv(i+1) yv(i)];
        zp = [zv(i) zv(i+1)       0       0];
        patch(xp, yp, zp, col, 'FaceAlpha', alph, 'EdgeColor', 'none')
    end
end

function C3_matrix = loadMatrix_fromData(data)
    C1_vals = data(:,1);
    C2_vals = data(:,2);
    C3_vals = data(:,3);
    C1_values = unique(C1_vals);
    C2_values = unique(C2_vals);
    C3_matrix = zeros(length(C1_values), length(C2_values));
    for k = 1:length(C3_vals)
        iC1 = find(C1_values == C1_vals(k));
        iC2 = find(C2_values == C2_vals(k));
        C3_matrix(iC1, iC2) = C3_vals(k);
    end
end

function data = loadAllFiles(files)
    data = [];
    for i = 1:length(files)
        opts = detectImportOptions(files{i}, 'FileType','text');
        opts.CommentStyle = '#'; opts.Delimiter = ' ';
        opts = setvartype(opts, 'double');
        data = [data; readmatrix(files{i}, opts)];
    end
end