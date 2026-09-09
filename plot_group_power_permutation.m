function plot_group_power_permutation(group_powA, group_powB, t_centers, cleanA, cleanB, state_name, chanlocs, n_perms, alpha_level, output_dir)
    
    num_subjs = length(group_powA);
    num_channels = size(group_powA{1}, 1);
    num_windows = length(t_centers);
    num_tests = num_channels * num_windows;
    
    subj_win_A = zeros(num_subjs, num_channels, num_windows);
    subj_win_B = zeros(num_subjs, num_channels, num_windows);
    trial_win_A = cell(num_subjs, 1);
    trial_win_B = cell(num_subjs, 1);
    
    valid_subjs = 0;
    
    % --- 1. Subject-Level Normalization ---
    for s = 1:num_subjs
        if isempty(group_powA{s}), continue; end
        valid_subjs = valid_subjs + 1;
        
        tmp_A = group_powA{s}; % [Chans x Windows x Trials]
        tmp_B = group_powB{s};
        
        % Z-score across pooled task windows to handle scale discrepancies
        subj_pool = cat(3, tmp_A, tmp_B);
        subj_mu   = mean(subj_pool(:), 'omitnan');
        subj_sd   = std(subj_pool(:), 0, 'omitnan');
        if subj_sd == 0 || isnan(subj_sd), subj_sd = 1; end
        
        tmp_A_norm = (tmp_A - subj_mu) / subj_sd;
        tmp_B_norm = (tmp_B - subj_mu) / subj_sd;
        
        trial_win_A{s} = tmp_A_norm;
        trial_win_B{s} = tmp_B_norm;
        
        subj_win_A(s, :, :) = mean(tmp_A_norm, 3, 'omitnan');
        subj_win_B(s, :, :) = mean(tmp_B_norm, 3, 'omitnan');
    end
    
    if valid_subjs < 2, return; end
    
    % --- 2. Compute Observed Group Differences ---
    group_A = squeeze(mean(subj_win_A, 1, 'omitnan')); % Row 1 
    group_B = squeeze(mean(subj_win_B, 1, 'omitnan')); % Row 2
    
    obs_diff = group_B - group_A; % Row 3
    obs_abs_diff = abs(obs_diff);
    
    % --- 3. Max-Statistic Permutation Testing (FWER) ---
    max_null_distribution = zeros(n_perms, 1);
    for p = 1:n_perms
        pseudo_subj_A = zeros(num_subjs, num_channels, num_windows);
        pseudo_subj_B = zeros(num_subjs, num_channels, num_windows);
        
        for s = 1:num_subjs
            if isempty(trial_win_A{s}), continue; end
            all_trials = cat(3, trial_win_A{s}, trial_win_B{s});
            n_total = size(all_trials, 3);
            nA = size(trial_win_A{s}, 3);
            
            shuffled_idx = randperm(n_total);
            
            pseudo_subj_A(s, :, :) = mean(all_trials(:, :, shuffled_idx(1:nA)), 3, 'omitnan');
            pseudo_subj_B(s, :, :) = mean(all_trials(:, :, shuffled_idx(nA+1:end)), 3, 'omitnan');
        end
        
        pseudo_grp_A = squeeze(mean(pseudo_subj_A, 1, 'omitnan'));
        pseudo_grp_B = squeeze(mean(pseudo_subj_B, 1, 'omitnan'));
        pseudo_diff  = pseudo_grp_B - pseudo_grp_A;
        
        max_null_distribution(p) = max(abs(pseudo_diff(:)));
    end
    
    % --- 4. Statistical Significance Thresholding & Printout ---
    threshold_val = prctile(max_null_distribution, (1 - alpha_level) * 100);
    sig_mask = obs_abs_diff > threshold_val;
    
    % Print significant channels to terminal
    total_sig = sum(sig_mask(:));
    if total_sig > 0
        fprintf('    -> FOUND %d SIGNIFICANT CHANNELS:\n', total_sig);
        for w = 1:num_windows
            if any(sig_mask(:, w))
                sig_ch_names = {chanlocs(sig_mask(:, w)).labels};
                fprintf('       [%5.2fs]: %s\n', t_centers(w), strjoin(sig_ch_names, ', '));
            end
        end
    else
        fprintf('    -> No significant channels found for %s.\n', state_name);
    end
    
    % --- 5. Plotting Topographies ---
    fig_width = max(1300, 220 * num_windows + 100);
    fig = figure('Position', [50, 80, fig_width, 860], 'Name', state_name, 'Visible', 'off');
    t = tiledlayout(3, num_windows, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    max_shared = max([abs(group_A(:)); abs(group_B(:)); abs(obs_diff(:))]);
    if max_shared == 0 || isnan(max_shared), max_shared = 0.3; end
    clim_Shared = [-max_shared, max_shared];
    
    row_titles = {sprintf('%s (Avg)', cleanA), sprintf('%s (Avg)', cleanB), sprintf('Diff (%s - %s)', cleanB, cleanA)};
    
    for row = 1:3
        for col = 1:num_windows
            nexttile;
            
            if row == 1
                data_to_plot = group_A(:, col);
            elseif row == 2
                data_to_plot = group_B(:, col);
            elseif row == 3
                data_to_plot = obs_diff(:, col);
            end
            
            if row == 3 && any(sig_mask(:, col))
                sig_chans = find(sig_mask(:, col));
                topoplot(data_to_plot, chanlocs, 'maplimits', clim_Shared, 'colormap', jet, ...
                    'emarker2', {sig_chans, '.', 'k', 16, 2});
            else
                topoplot(data_to_plot, chanlocs, 'maplimits', clim_Shared, 'colormap', jet);
            end
            
            if col == 1
                text(-0.72, 0, row_titles{row}, 'FontSize', 16, 'FontWeight', 'bold', ...
                    'Rotation', 90, 'HorizontalAlignment', 'center');
            end
            
            if row == 3
                text(0, -0.68, sprintf('%.2fs', t_centers(col)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
            end
            
            % NEW: Attach colorbars natively to the last tile of each row to prevent overlapping
            if col == num_windows
                cb = colorbar;
                cb.FontSize = 14;
                if row == 3
                    ylabel(cb, '\Delta z', 'FontSize', 16, 'FontWeight', 'bold');
                else
                    ylabel(cb, 'Normalized z', 'FontSize', 16, 'FontWeight', 'bold');
                end
            end
        end
    end
    
    sgtitle(sprintf('Group-Level Permutation Test: %s', state_name), 'FontSize', 22, 'FontWeight', 'bold');
    
    % Safely strip out spaces, parentheses, AND forward slashes for the filename
    safe_name = strrep(strrep(strrep(strrep(state_name, ' ', '_'), '(', ''), ')', ''), '/', '_');
    
    save_name = fullfile(output_dir, sprintf('Group_Power_Perm_%s_vs_%s_%s.png', strrep(cleanA,' ',''), strrep(cleanB,' ',''), safe_name));
    saveas(fig, save_name);
    close(fig);
end