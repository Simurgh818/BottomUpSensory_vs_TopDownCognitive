function plot_dynamic_networks(sub_ch_A_all, sub_ch_B_all, sub_pc_A_all, sub_pc_B_all, window_centers, ...
    trialsA, trialsB, time_ms_eeg, t_win, fs, window_size_ms, step_size_ms, ...
    condA, condB, active_state, current_band, k_opt, pc_labels, all_channels_str, output_dir, subj_id)
    
    num_windows = length(window_centers);
    if num_windows == 0, return; end
    num_ch = size(sub_ch_A_all, 1);
    
    get_clean = @(c) strrep(strrep(strrep(strrep(strrep(strrep(c, ...
        'BLA', 'Auditory'), 'BLT', 'Tactile'), 'P1', 'Cued'), ...
        'P2', 'Unpred.'), 'P3', 'Rand. Cued'), '_', ' ');
    cleanA = get_clean(condA);
    cleanB = get_clean(condB);
    
    % --- RESTORED: Dynamic Y-Axis Limits (Channels and PCs) ---
    switch lower(current_band)
        case 'delta', ch_ylim = [-2.5, 2.5]; pc_ylim = [-5, 5];
        case 'theta', ch_ylim = [-2.6, 2.6]; pc_ylim = [-7.5, 7.5];
        case 'alpha', ch_ylim = [-1.5, 1.5]; pc_ylim = [-3, 3];
        case 'beta',  ch_ylim = [-1, 1];     pc_ylim = [-2, 2];
        case 'gamma', ch_ylim = [-1, 1];     pc_ylim = [-2, 2];
        otherwise,    ch_ylim = [-5, 5];     pc_ylim = [-15, 15];
    end
    if strcmp(condA, 'BLA') && strcmp(condB, 'P1')
        ch_ylim = ch_ylim * 2;
        pc_ylim = pc_ylim * 2;
    end
    
    idx_start = find(time_ms_eeg >= t_win(1), 1, 'first');
    win_samples = round((window_size_ms / 1000) * fs);
    step_samples = round((step_size_ms / 1000) * fs);
    start_idx = 1 : step_samples : (length(time_ms_eeg(idx_start:find(time_ms_eeg <= t_win(2), 1, 'last'))) - win_samples + 1);

    corr_thresh = 0.5; 
    c_lim_mat   = [-1 1]; 
    cmap = jet(256);
    get_color = @(v) cmap(max(1, min(256, round((v + 1) / 2 * 255) + 1)), :);
    get_high_corr_idx = @(mat, thresh) find(any(abs(mat - diag(diag(mat))) > thresh, 2));

    %% =========================================================================
    % FIGURE 1: MASTER CHANNELS (Heatmaps, Networks & Traces)
    % =========================================================================
    active_channels = [];
    for w = 1:num_windows
        active_channels = union(active_channels, union(get_high_corr_idx(sub_ch_A_all(:,:,w), corr_thresh), ...
                                                       get_high_corr_idx(sub_ch_B_all(:,:,w), corr_thresh)));
    end
    ch_colors = hsv(num_ch); 
    
    if exist('EEG', 'var') && isfield(EEG, 'chanlocs')
        th = pi/180 * [EEG.chanlocs.theta]; rd = [EEG.chanlocs.radius];
        x_node = rd .* sin(th); y_node = rd .* cos(th); rmax = max(0.5, max(rd)); 
    else
        th = linspace(0, 2*pi, num_ch+1); th(end)=[];
        x_node = 0.5 * cos(th); y_node = 0.5 * sin(th); rmax = 0.5;
    end
    
    fig_width = max(1200, 200 * num_windows);
    fig1 = figure('Position', [50, 50, fig_width, 1600], 'Name', sprintf('[%s] %s vs %s Master Network', upper(current_band), cleanA, cleanB), 'Visible', 'off');
    tiledlayout(6, num_windows, 'TileSpacing', 'tight', 'Padding', 'tight');
    
    theta_circle = linspace(0, 2*pi, 100);
    draw_head = @() plot(rmax*cos(theta_circle), rmax*sin(theta_circle), 'k', 'LineWidth', 1.5);
    draw_nose = @() plot([rmax*0.1, 0, -rmax*0.1], [rmax, rmax*1.15, rmax], 'k', 'LineWidth', 1.5);
    
    % ROW 1: Cond A Heatmap
    for w = 1:num_windows
        nexttile; imagesc(sub_ch_A_all(:,:,w)); colormap(gca, 'jet'); clim(c_lim_mat); axis square;
        xticks(1:num_ch); xticklabels({}); 
        if w == 1
            ylabel(sprintf('%s\n(Channels)', cleanA), 'FontSize', 16, 'FontWeight', 'bold');
            yticks(1:num_ch); yticklabels(all_channels_str); ax = gca; ax.YAxis.FontSize = 8;
        else yticks([]); end
        title(sprintf('%.2fs', window_centers(w)), 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    % ROW 2: Cond B Heatmap
    for w = 1:num_windows
        nexttile; imagesc(sub_ch_B_all(:,:,w)); colormap(gca, 'jet'); clim(c_lim_mat); axis square;
        xticks(1:num_ch); xticklabels(all_channels_str); xtickangle(90); ax = gca; ax.XAxis.FontSize = 8;
        if w == 1
            ylabel(sprintf('%s\n(Channels)', cleanB), 'FontSize', 16, 'FontWeight', 'bold');
            yticks(1:num_ch); yticklabels(all_channels_str); ax = gca; ax.YAxis.FontSize = 8;
        else yticks([]); end
        if w == num_windows, ax_master_color = gca; end
    end
    
    % ROW 3: Cond A Network Graph
    for w = 1:num_windows
        nexttile; hold on; draw_head(); draw_nose(); scatter(x_node, y_node, 15, 'k', 'filled'); 
        matA = sub_ch_A_all(:,:,w); matA(1:num_ch+1:end) = 0; [row, col] = find(abs(matA) > corr_thresh);
        for i = 1:length(row)
            if row(i) > col(i) 
                v = matA(row(i), col(i));
                plot([x_node(row(i)) x_node(col(i))], [y_node(row(i)) y_node(col(i))], 'Color', get_color(v), 'LineWidth', abs(v) * 4);
            end
        end
        axis equal; axis off; colormap(gca, 'jet'); clim(c_lim_mat); 
        if w == 1, text(-rmax*1.5, 0, sprintf('%s\n(|\\Delta r| > %g)', cleanA, corr_thresh), 'Rotation', 90, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 16); end
    end
    
    % ROW 4: Cond B Network Graph
    for w = 1:num_windows
        nexttile; hold on; draw_head(); draw_nose(); scatter(x_node, y_node, 15, 'k', 'filled'); 
        matB = sub_ch_B_all(:,:,w); matB(1:num_ch+1:end) = 0; [row, col] = find(abs(matB) > corr_thresh);
        for i = 1:length(row)
            if row(i) > col(i) 
                v = matB(row(i), col(i));
                plot([x_node(row(i)) x_node(col(i))], [y_node(row(i)) y_node(col(i))], 'Color', get_color(v), 'LineWidth', abs(v) * 4);
            end
        end
        axis equal; axis off; colormap(gca, 'jet'); clim(c_lim_mat);
        if w == 1, text(-rmax*1.5, 0, sprintf('%s\n(|\\Delta r| > %g)', cleanB, corr_thresh), 'Rotation', 90, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 16); end
    end

    % ROW 5: Cond A Traces
    for w = 1:num_windows
        nexttile; hold on; w_s = idx_start + start_idx(w) - 1; w_e = w_s + win_samples - 1; t_axis = time_ms_eeg(w_s:w_e);
        idx = get_high_corr_idx(sub_ch_A_all(:,:,w), corr_thresh);
        if ~isempty(idx)
            for i = 1:length(idx), plot(t_axis, mean(trialsA(idx(i), w_s:w_e, :), 3, 'omitnan'), 'LineWidth', 1.2, 'Color', ch_colors(idx(i), :)); end
        else
            text(mean(t_axis), 0, sprintf('None (> %g)', corr_thresh), 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
        end
        xlim([t_axis(1), t_axis(end)]); ylim(ch_ylim); set(gca, 'XColor', 'none'); 
        if w == 1, ylabel(sprintf('%s\nCh Traces', cleanA), 'FontWeight', 'bold', 'FontSize', 16); end
        
        if w == num_windows && ~isempty(active_channels)
            h_leg = arrayfun(@(x) plot(nan, nan, 'LineWidth', 1.5, 'Color', ch_colors(x, :)), active_channels);
            lgd = legend(h_leg, all_channels_str(active_channels), 'Location', 'eastoutside', 'FontSize', 8, 'AutoUpdate', 'off'); lgd.ItemTokenSize = [12, 12];
        end
    end
    
    % ROW 6: Cond B Traces
    for w = 1:num_windows
        nexttile; hold on; w_s = idx_start + start_idx(w) - 1; w_e = w_s + win_samples - 1; t_axis = time_ms_eeg(w_s:w_e);
        idx = get_high_corr_idx(sub_ch_B_all(:,:,w), corr_thresh);
        if ~isempty(idx)
            for i = 1:length(idx), plot(t_axis, mean(trialsB(idx(i), w_s:w_e, :), 3, 'omitnan'), 'LineWidth', 1.2, 'Color', ch_colors(idx(i), :)); end
        else
            text(mean(t_axis), 0, sprintf('None (> %g)', corr_thresh), 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
        end
        xlim([t_axis(1), t_axis(end)]); ylim(ch_ylim); 
        if w == 1, ylabel(sprintf('%s\nCh Traces', cleanB), 'FontWeight', 'bold', 'FontSize', 16); end
        xlabel('Time', 'FontSize', 16, 'FontWeight', 'bold');
    end
    
    cb = colorbar(ax_master_color); cb.Layout.Tile = 'east'; cb.Label.String = 'Subtracted Correlation (\Delta r)'; cb.Label.FontSize = 16; cb.Label.FontWeight = 'bold';
    sgtitle(sprintf('[%s Band] Channels: Heatmaps, Node-Link & Traces (%s vs %s)', upper(current_band), cleanA, cleanB), 'FontSize', 22, 'FontWeight', 'bold');
    saveas(fig1, fullfile(output_dir, sprintf('%s_%s_vs_%s_Channels_%s.png', subj_id, condA, condB, strrep(active_state,' ','_'))));
    close(fig1);

    %% =========================================================================
    % FIGURE 2: MASTER LATENT dPCs (Heatmaps, Trajectories & Traces)
    % =========================================================================
    % 1. Local SVD projection into Shared Latent Space
    avgA_raw = mean(trialsA, 3, 'omitnan'); avgB_raw = mean(trialsB, 3, 'omitnan');
    dat_comb = [avgA_raw, avgB_raw] - mean([avgA_raw, avgB_raw], 2);
    [~, ~, V] = svd(dat_comb', 'econ'); W_proj = V(:, 1:k_opt)';
    
    datA_pc = reshape(W_proj * reshape(trialsA, num_ch, []), k_opt, size(trialsA,2), size(trialsA,3));
    avgA_pc = mean(datA_pc, 3, 'omitnan');
    
    datB_pc = reshape(W_proj * reshape(trialsB, num_ch, []), k_opt, size(trialsB,2), size(trialsB,3));
    avgB_pc = mean(datB_pc, 3, 'omitnan');

    active_pcs = [];
    for w = 1:num_windows
        active_pcs = union(active_pcs, union(get_high_corr_idx(sub_pc_A_all(:,:,w), corr_thresh), ...
                                             get_high_corr_idx(sub_pc_B_all(:,:,w), corr_thresh)));
    end
    pc_colors = lines(k_opt); 

    fig2 = figure('Position', [100, 100, fig_width, 1600], 'Name', sprintf('[%s] %s vs %s Master dPCs', upper(current_band), cleanA, cleanB), 'Visible', 'off');
    tiledlayout(6, num_windows, 'TileSpacing', 'tight', 'Padding', 'tight');
    
    % ROW 1: Cond A dPC Heatmap
    for w = 1:num_windows
        nexttile; imagesc(sub_pc_A_all(:,:,w)); colormap(gca, 'jet'); clim(c_lim_mat); axis square;
        xticks(1:k_opt); xticklabels({}); 
        if w == 1
            ylabel(sprintf('%s\n(Shared dPCs)', cleanA), 'FontSize', 16, 'FontWeight', 'bold');
            yticks(1:k_opt); yticklabels(pc_labels); ax = gca; ax.YAxis.FontSize = 8;
        else yticks([]); end
        title(sprintf('%.2fs', window_centers(w)), 'FontSize', 14, 'FontWeight', 'bold');
    end

    % ROW 2: Cond B dPC Heatmap
    for w = 1:num_windows
        nexttile; imagesc(sub_pc_B_all(:,:,w)); colormap(gca, 'jet'); clim(c_lim_mat); axis square;
        xticks(1:k_opt); xticklabels(pc_labels); xtickangle(45); ax = gca; ax.XAxis.FontSize = 8;
        if w == 1
            ylabel(sprintf('%s\n(Shared dPCs)', cleanB), 'FontSize', 16, 'FontWeight', 'bold');
            yticks(1:k_opt); yticklabels(pc_labels); ax = gca; ax.YAxis.FontSize = 8;
        else yticks([]); end
        if w == num_windows, ax_dpc_color = gca; end
    end

    % --- ROW 3: Cond A State-Space Trajectory (dPC1 vs dPC2) ---
    for w = 1:num_windows
        nexttile; hold on;
        w_s = idx_start + start_idx(w) - 1; w_e = w_s + win_samples - 1;
        plot(avgA_pc(1, w_s:w_e), avgA_pc(2, w_s:w_e), 'LineWidth', 2, 'Color', [0.2 0.6 0.8]);
        scatter(avgA_pc(1, w_s), avgA_pc(2, w_s), 40, [0.2 0.6 0.8], 'filled'); % Dot indicates start of time window
        xlim(pc_ylim); ylim(pc_ylim); grid on;
        
        if w == 1
            % UPDATED: Y-Label changed to clearly state 'dPC 2'
            ylabel(sprintf('%s\ndPC 2', cleanA), 'FontSize', 16, 'FontWeight', 'bold');
        else 
            yticks([]); 
        end
        xticks([]); % Hide X-ticks for Row 3 to keep it flush with Row 4
    end

    % --- ROW 4: Cond B State-Space Trajectory (dPC1 vs dPC2) ---
    for w = 1:num_windows
        nexttile; hold on;
        w_s = idx_start + start_idx(w) - 1; w_e = w_s + win_samples - 1;
        plot(avgB_pc(1, w_s:w_e), avgB_pc(2, w_s:w_e), 'LineWidth', 2, 'Color', [0.8 0.4 0.2]);
        scatter(avgB_pc(1, w_s), avgB_pc(2, w_s), 40, [0.8 0.4 0.2], 'filled'); % Dot indicates start of time window
        xlim(pc_ylim); ylim(pc_ylim); grid on;
        
        if w == 1
            % UPDATED: Y-Label changed to clearly state 'dPC 2'
            ylabel(sprintf('%s\ndPC 2', cleanB), 'FontSize', 16, 'FontWeight', 'bold');
        else 
            yticks([]); 
        end
        
        % UPDATED: Removed xticks([]) to reveal the numerical scale, and added the X-Label
        xlabel('dPC 1', 'FontSize', 16, 'FontWeight', 'bold');
    end

    % ROW 5: Cond A dPC Traces (> 0.5 Corr)
    for w = 1:num_windows
        nexttile; hold on; w_s = idx_start + start_idx(w) - 1; w_e = w_s + win_samples - 1; t_axis = time_ms_eeg(w_s:w_e);
        idx = get_high_corr_idx(sub_pc_A_all(:,:,w), corr_thresh);
        if ~isempty(idx)
            for i = 1:length(idx), plot(t_axis, avgA_pc(idx(i), w_s:w_e), 'LineWidth', 1.2, 'Color', pc_colors(idx(i), :)); end
        else
            text(mean(t_axis), 0, sprintf('None (> %g)', corr_thresh), 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
        end
        xlim([t_axis(1), t_axis(end)]); ylim(pc_ylim); set(gca, 'XColor', 'none'); 
        if w == 1, ylabel(sprintf('%s\ndPC Traces', cleanA), 'FontWeight', 'bold', 'FontSize', 16); end
        
        if w == num_windows && ~isempty(active_pcs)
            h_leg = arrayfun(@(x) plot(nan, nan, 'LineWidth', 1.5, 'Color', pc_colors(x, :)), active_pcs);
            lgd = legend(h_leg, pc_labels(active_pcs), 'Location', 'eastoutside', 'FontSize', 8, 'AutoUpdate', 'off'); lgd.ItemTokenSize = [12, 12];
        end
    end

    % ROW 6: Cond B dPC Traces (> 0.5 Corr)
    for w = 1:num_windows
        nexttile; hold on; w_s = idx_start + start_idx(w) - 1; w_e = w_s + win_samples - 1; t_axis = time_ms_eeg(w_s:w_e);
        idx = get_high_corr_idx(sub_pc_B_all(:,:,w), corr_thresh);
        if ~isempty(idx)
            for i = 1:length(idx), plot(t_axis, avgB_pc(idx(i), w_s:w_e), 'LineWidth', 1.2, 'Color', pc_colors(idx(i), :)); end
        else
            text(mean(t_axis), 0, sprintf('None (> %g)', corr_thresh), 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
        end
        xlim([t_axis(1), t_axis(end)]); ylim(pc_ylim); 
        if w == 1, ylabel(sprintf('%s\ndPC Traces', cleanB), 'FontWeight', 'bold', 'FontSize', 16); end
        xlabel('Time', 'FontSize', 16, 'FontWeight', 'bold');
    end

    cb2 = colorbar(ax_dpc_color); cb2.Layout.Tile = 'east'; cb2.Label.String = 'Subtracted Correlation (\Delta r)'; cb2.Label.FontSize = 16; cb2.Label.FontWeight = 'bold';
    sgtitle(sprintf('[%s Band] Latent Subspace: Heatmaps, State-Space & Traces (%s vs %s)', upper(current_band), cleanA, cleanB), 'FontSize', 22, 'FontWeight', 'bold');
    
    saveas(fig2, fullfile(output_dir, sprintf('%s_%s_vs_%s_dPCs_%s.png', subj_id, condA, condB, strrep(active_state,' ','_'))));
    close(fig2);
end