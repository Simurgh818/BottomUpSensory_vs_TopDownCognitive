function plot_correlation_filmstrip(trials_filt, t_win, time_ms_eeg, window_size_ms, step_size_ms, bg_windows, clean_cond, band_name, all_channels_str, output_dir, subj_id)
    num_ch = size(trials_filt, 1);
    num_trials = size(trials_filt, 3);
    if num_trials < 2, return; end
    
    % --- 1. Bootstrapped Resting-State Baseline & Null Distribution ---
    idx_bg = find(time_ms_eeg >= bg_windows(1) - 1e-5 & time_ms_eeg <= bg_windows(2) + 1e-5);
    trials_bg = trials_filt(:, idx_bg, :);
    
    n_boot = 1000;
    bg_dist = zeros(num_ch, num_ch, n_boot);
    
    for i = 1:n_boot
        samp_idx = randi(num_trials, num_trials, 1);
        avg_boot = mean(trials_bg(:, :, samp_idx), 3, 'omitnan');
        bg_dist(:,:,i) = corrcoef(avg_boot');
    end
    
    bg_ch_stable = mean(bg_dist, 3, 'omitnan');
    delta_bg_dist = abs(bg_dist - bg_ch_stable);
    threshold_matrix = prctile(delta_bg_dist, 95, 3);
    
    % --- 2. Sliding Window Active Correlation ---
    window_size_s = window_size_ms / 1000;
    step_size_s = step_size_ms / 1000;
    t_centers = (t_win(1) + step_size_s/2) : step_size_s : (t_win(2) - step_size_s/2 + 1e-5);
    num_windows = length(t_centers);
    
    avg_act = mean(trials_filt, 3, 'omitnan');
    raw_ch_all = nan(num_ch, num_ch, num_windows);
    sub_ch_all = nan(num_ch, num_ch, num_windows);
    
    for w = 1:num_windows
        wc = t_centers(w);
        idx_act = find(time_ms_eeg >= (wc - window_size_s/2 - 1e-5) & time_ms_eeg <= (wc + window_size_s/2 + 1e-5));
        
        real_ch = corrcoef(avg_act(:, idx_act)');
        sub_ch  = real_ch - bg_ch_stable;
        
        % 95% Null Distribution Masking
        sig_mask   = abs(sub_ch) > threshold_matrix;
        sub_ch_sig = sub_ch .* sig_mask;
        
        raw_ch_all(:,:,w) = real_ch;
        sub_ch_all(:,:,w) = sub_ch_sig;
    end
    
    % --- 3. 3-Row Filmstrip Plotting ---
    c_lim = [-1 1];
    fig_width = max(1000, 220 * num_windows);
    fig = figure('Position', [50, 50, fig_width, 750], 'Name', sprintf('%s %s Filmstrip', clean_cond, band_name), 'Visible', 'off');
    tiledlayout(3, num_windows, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    % Row 1: Raw Active Correlation
    for w = 1:num_windows
        nexttile; imagesc(raw_ch_all(:,:,w));
        colormap('jet'); clim(c_lim); axis square; set(gca, 'FontSize', 8);
        if w == 1
            text(-0.35, 0.5, sprintf('%s (%s)\nActive Window', clean_cond, upper(band_name)), 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold');
            yticks(1:num_ch); yticklabels(all_channels_str);
        else
            yticks([]);
        end
        xticks([]);
        w_start = t_centers(w) - (step_size_s / 2);
        title(sprintf('%.2fs', w_start), 'FontSize', 15, 'FontWeight', 'bold');
    end
    
    % Row 2: Bootstrapped Baseline Matrix
    for w = 1:num_windows
        nexttile; imagesc(bg_ch_stable);
        colormap('jet'); clim(c_lim); axis square; set(gca, 'FontSize', 8);
        if w == 1
            text(-0.35, 0.5, sprintf('Baseline (0.0-0.5s)\n1000 Bootstraps'), 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold');
            yticks(1:num_ch); yticklabels(all_channels_str);
        else
            yticks([]);
        end
        xticks([]);
    end
    
    % Row 3: Subtracted & 95% Thresholded Matrix
    for w = 1:num_windows
        nexttile; imagesc(sub_ch_all(:,:,w));
        colormap('jet'); clim(c_lim); axis square; set(gca, 'FontSize', 8);
        if w == 1
            text(-0.35, 0.5, sprintf('Sig. \\Delta r\n(>95%% Null)'), 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold');
            yticks(1:num_ch); yticklabels(all_channels_str);
        else
            yticks([]);
        end
        xticks(1:num_ch); xticklabels(all_channels_str); xtickangle(90);
    end
    
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Label.String = 'Correlation (r)';
    cb.Label.FontSize = 14;
    cb.Label.FontWeight = 'bold';
    
    sgtitle(sprintf('%s [%s Band]: Raw, Baseline, and Significant Subtracted Networks', clean_cond, upper(band_name)), ...
        'FontSize', 20, 'FontWeight', 'bold');
    
    save_file = fullfile(output_dir, sprintf('%s_%s_%s_CorrFilmstrip.png', subj_id, strrep(clean_cond,' ',''), upper(band_name)));
    saveas(fig, save_file);
    close(fig);
end