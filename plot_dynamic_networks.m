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

    c_lim_mat = [-1 1]; % Using strictly normalized correlation [-1, 1]
    fig_width = max(1000, 250 * num_windows);
    
    %% --- FIGURE 1: HEATMAPS ---
    fig1 = figure('Position', [50, 50, fig_width, 1000], 'Visible', 'off');
    tiledlayout(4, num_windows, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for w = 1:num_windows, nexttile; imagesc(sub_ch_A_all(:,:,w)); colormap('jet'); clim(c_lim_mat); axis square; xticks([]); if w==1, ylabel(sprintf('%s\n(Channels)', cleanA), 'FontWeight', 'bold'); yticks(1:num_ch); yticklabels(all_channels_str); else yticks([]); end; title(sprintf('%.2fs', window_centers(w))); end
    for w = 1:num_windows, nexttile; imagesc(sub_ch_B_all(:,:,w)); colormap('jet'); clim(c_lim_mat); axis square; xticks([]); if w==1, ylabel(sprintf('%s\n(Channels)', cleanB), 'FontWeight', 'bold'); yticks(1:num_ch); yticklabels(all_channels_str); else yticks([]); end; end
    for w = 1:num_windows, nexttile; imagesc(sub_pc_A_all(:,:,w)); colormap('jet'); clim(c_lim_mat); axis square; xticks([]); if w==1, ylabel(sprintf('%s\n(Shared dPCs)', cleanA), 'FontWeight', 'bold'); yticks(1:k_opt); yticklabels(pc_labels); else yticks([]); end; end
    for w = 1:num_windows, nexttile; imagesc(sub_pc_B_all(:,:,w)); colormap('jet'); clim(c_lim_mat); axis square; xticks(1:k_opt); xticklabels(pc_labels); xtickangle(45); if w==1, ylabel(sprintf('%s\n(Shared dPCs)', cleanB), 'FontWeight', 'bold'); yticks(1:k_opt); yticklabels(pc_labels); else yticks([]); end; end
    
    cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'Subtracted Correlation (\Delta r)'; cb.Label.FontWeight = 'bold';
    sgtitle(sprintf('[%s Band] Network Evolution (%s): %s vs %s', upper(current_band), active_state, cleanA, cleanB), 'FontSize', 18, 'FontWeight', 'bold');
    
    % DYNAMIC SAVE: e.g., BOS10_BLA_vs_P1_Heatmaps_Stimulus_State.png
    saveas(fig1, fullfile(output_dir, sprintf('%s_%s_vs_%s_Heatmaps_%s.png', subj_id, condA, condB, strrep(active_state,' ','_'))));
    close(fig1);

    %% --- FIGURE 2: FILTERED TRACES ---
    corr_thresh = 0.5;
    fig2 = figure('Position', [100, 100, fig_width, 900], 'Visible', 'off');
    tiledlayout(4, num_windows, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Dynamic limits
    if strcmp(condA, 'BLA') && strcmp(condB, 'P1')
        ch_ylim = [-10, 10]; pc_ylim = [-25, 25];
    else
        ch_ylim = [-5, 5]; pc_ylim = [-15, 15];
    end
    
    idx_start = find(time_ms_eeg >= t_win(1), 1, 'first');
    win_samples = round((window_size_ms / 1000) * fs);
    step_samples = round((step_size_ms / 1000) * fs);
    start_idx = 1 : step_samples : (length(time_ms_eeg(idx_start:find(time_ms_eeg <= t_win(2), 1, 'last'))) - win_samples + 1);

    get_idx = @(mat, thresh) find(any(abs(mat - diag(diag(mat))) > thresh, 2));
    ch_c = hsv(num_ch); pc_c = lines(k_opt);
    
    active_ch = []; active_pc = [];
    for w = 1:num_windows
        active_ch = union(active_ch, union(get_idx(sub_ch_A_all(:,:,w), corr_thresh), get_idx(sub_ch_B_all(:,:,w), corr_thresh)));
        active_pc = union(active_pc, union(get_idx(sub_pc_A_all(:,:,w), corr_thresh), get_idx(sub_pc_B_all(:,:,w), corr_thresh)));
    end
    
    for row = 1:4
        for w = 1:num_windows
            nexttile; hold on;
            w_s = idx_start + start_idx(w) - 1; w_e = w_s + win_samples - 1; t_axis = time_ms_eeg(w_s:w_e);
            
            if row == 1, mat = sub_ch_A_all(:,:,w); dat = trialsA; col = ch_c; yl = ch_ylim;
            elseif row == 2, mat = sub_ch_B_all(:,:,w); dat = trialsB; col = ch_c; yl = ch_ylim;
            elseif row == 3, mat = sub_pc_A_all(:,:,w); dat = trialsA; col = pc_c; yl = pc_ylim; 
            else, mat = sub_pc_B_all(:,:,w); dat = trialsB; col = pc_c; yl = pc_ylim; end
            
            % Local projection for plotting dPCs
            if row >= 3
                avgA_raw = mean(trialsA, 3, 'omitnan'); avgB_raw = mean(trialsB, 3, 'omitnan');
                dat_comb = [avgA_raw, avgB_raw] - mean([avgA_raw, avgB_raw], 2);
                [~, ~, V] = svd(dat_comb', 'econ'); W_proj = V(:, 1:k_opt)';
                dat_flat = reshape(dat, num_ch, []);
                dat = reshape(W_proj * dat_flat, k_opt, size(dat,2), size(dat,3));
            end
            
            idx = get_idx(mat, corr_thresh);
            if ~isempty(idx)
                for i=1:length(idx), plot(t_axis, mean(dat(idx(i), w_s:w_e, :), 3), 'Color', col(idx(i),:), 'LineWidth', 1.2); end
            else
                text(mean(t_axis), 0, sprintf('None > %g', corr_thresh), 'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'center');
            end
            xlim([t_axis(1) t_axis(end)]); ylim(yl); set(gca, 'XColor', 'none');
            
            if w == 1
                if row == 1, ylabel(sprintf('%s\nCh Traces', cleanA), 'FontWeight', 'bold'); title(sprintf('%.2fs', window_centers(w)));
                elseif row == 2, ylabel(sprintf('%s\nCh Traces', cleanB), 'FontWeight', 'bold');
                elseif row == 3, ylabel(sprintf('%s\ndPC Traces', cleanA), 'FontWeight', 'bold');
                else, ylabel(sprintf('%s\ndPC Traces', cleanB), 'FontWeight', 'bold'); xlabel('Time'); end
            elseif row == 1
                title(sprintf('%.2fs', window_centers(w)), 'FontSize', 12, 'FontWeight', 'bold');
            end
            
            % Add Unified Legends to rows 2 & 4
            if w == num_windows
                if row == 2 && ~isempty(active_ch)
                    h_leg = arrayfun(@(x) plot(nan, nan, 'Color', ch_c(x,:), 'LineWidth', 1.5), active_ch);
                    lgd = legend(h_leg, all_channels_str(active_ch), 'Location', 'eastoutside', 'FontSize', 8, 'AutoUpdate', 'off'); lgd.ItemTokenSize = [12, 12];
                elseif row == 4 && ~isempty(active_pc)
                    h_leg = arrayfun(@(x) plot(nan, nan, 'Color', pc_c(x,:), 'LineWidth', 1.5), active_pc);
                    lgd = legend(h_leg, pc_labels(active_pc), 'Location', 'eastoutside', 'FontSize', 8, 'AutoUpdate', 'off'); lgd.ItemTokenSize = [12, 12];
                end
            end
        end
    end
    
    sgtitle(sprintf('[%s Band] Traces (%s): %s vs %s (>%g Corr)', upper(current_band), active_state, cleanA, cleanB, corr_thresh), 'FontSize', 18, 'FontWeight', 'bold');
    
    % DYNAMIC SAVE: e.g., BOS10_BLA_vs_P1_Traces_Stimulus_State.png
    saveas(fig2, fullfile(output_dir, sprintf('%s_%s_vs_%s_Traces_%s.png', subj_id, condA, condB, strrep(active_state,' ','_'))));
    close(fig2);
end