function plot_group_power_permutation(pow_A, pow_B, t_axis, condA_name, condB_name, state_name, chanlocs, n_perms, alpha_level, output_dir)
    % PLOT_GROUP_POWER_PERMUTATION Performs group-level max-statistic permutation
    % testing on band power topographies and generates publication-ready figures.

    if nargin < 8 || isempty(n_perms), n_perms = 1000; end
    if nargin < 9 || isempty(alpha_level), alpha_level = 0.05; end
    if nargin < 10 || isempty(output_dir), output_dir = pwd; end
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    num_windows = length(t_axis);
    num_subjs = length(pow_A);
    
    % Determine channel count from first subject
    if iscell(pow_A)
        num_channels = size(pow_A{1}, 1);
    else
        num_channels = size(pow_A, 2);
    end
    num_tests = num_channels * num_windows;

    % --- 1. Channel Label Extraction from Passed chanlocs ---
    if ~isempty(chanlocs) && isstruct(chanlocs) && isfield(chanlocs, 'labels')
        chan_labels = {chanlocs.labels};
    else
        chan_labels = arrayfun(@(x) sprintf('Ch%d', x), 1:num_channels, 'UniformOutput', false);
    end

    % --- 2. Extract Data & Apply Subject-Level Normalization ---
    subj_win_A = zeros(num_subjs, num_channels, num_windows);
    subj_win_B = zeros(num_subjs, num_channels, num_windows);
    trial_win_A = cell(num_subjs, 1);
    trial_win_B = cell(num_subjs, 1);

    is_trial_level = iscell(pow_A) && (ndims(pow_A{1}) == 3) && (size(pow_A{1}, 3) > 1);

    for s = 1:num_subjs
        if iscell(pow_A)
            datA = pow_A{s};
            datB = pow_B{s};
        else
            datA = squeeze(pow_A(s, :, :));
            datB = squeeze(pow_B(s, :, :));
        end

        % Subject z-score standardization across pooled conditions
        subj_pool = cat(ndims(datA), datA, datB);
        subj_mu   = mean(subj_pool(:), 'omitnan');
        subj_sd   = std(subj_pool(:), 0, 'omitnan');
        if subj_sd == 0 || isnan(subj_sd), subj_sd = 1; end

        normA = (datA - subj_mu) / subj_sd;
        normB = (datB - subj_mu) / subj_sd;

        trial_win_A{s} = normA;
        trial_win_B{s} = normB;

        if is_trial_level
            subj_win_A(s, :, :) = mean(normA, 3, 'omitnan');
            subj_win_B(s, :, :) = mean(normB, 3, 'omitnan');
        else
            subj_win_A(s, :, :) = normA;
            subj_win_B(s, :, :) = normB;
        end
    end

    % --- 3. Observed Group Averages & Differences ---
    group_A = squeeze(mean(subj_win_A, 1, 'omitnan')); % [Chans x Windows]
    group_B = squeeze(mean(subj_win_B, 1, 'omitnan'));
    
    obs_diff = group_B - group_A;
    obs_abs_diff = abs(obs_diff);

    % --- 4. Max-Statistic Permutation Testing (FWER Controlled) ---
    max_null_distribution = zeros(n_perms, 1);

    for p = 1:n_perms
        pseudo_subj_A = zeros(num_subjs, num_channels, num_windows);
        pseudo_subj_B = zeros(num_subjs, num_channels, num_windows);

        for s = 1:num_subjs
            if is_trial_level
                all_trials = cat(3, trial_win_A{s}, trial_win_B{s});
                n_total = size(all_trials, 3);
                nA = size(trial_win_A{s}, 3);
                shuffled_idx = randperm(n_total);

                pseudo_subj_A(s, :, :) = mean(all_trials(:, :, shuffled_idx(1:nA)), 3, 'omitnan');
                pseudo_subj_B(s, :, :) = mean(all_trials(:, :, shuffled_idx(nA+1:end)), 3, 'omitnan');
            else
                % Paired condition-label coin flip per subject
                if rand() > 0.5
                    pseudo_subj_A(s, :, :) = subj_win_A(s, :, :);
                    pseudo_subj_B(s, :, :) = subj_win_B(s, :, :);
                else
                    pseudo_subj_A(s, :, :) = subj_win_B(s, :, :);
                    pseudo_subj_B(s, :, :) = subj_win_A(s, :, :);
                end
            end
        end

        pseudo_grp_A = squeeze(mean(pseudo_subj_A, 1, 'omitnan'));
        pseudo_grp_B = squeeze(mean(pseudo_subj_B, 1, 'omitnan'));
        pseudo_diff  = pseudo_grp_B - pseudo_grp_A;

        max_null_distribution(p) = max(abs(pseudo_diff(:)));
    end

    % --- 5. Significance Thresholding ---
    threshold_val = prctile(max_null_distribution, (1 - alpha_level) * 100);
    sig_mask = obs_abs_diff > threshold_val;

    % --- 6. Plotting Topographies ---
    fig_width = max(1300, 220 * num_windows + 100);
    fig = figure('Position', [50, 80, fig_width, 860], 'Name', state_name, 'Visible', 'off');
    t = tiledlayout(3, num_windows, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Shared symmetric color boundaries
    max_shared = max([abs(group_A(:)); abs(group_B(:)); abs(obs_diff(:))]);
    if max_shared == 0 || isnan(max_shared), max_shared = 0.5; end
    clim_Shared = [-max_shared, max_shared];

    row_titles = {sprintf('%s (Avg)', condA_name), ...
                  sprintf('%s (Avg)', condB_name), ...
                  sprintf('Diff (%s - %s)', condB_name, condA_name)};

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

            % Plot topomap using passed chanlocs
            if row == 3 && any(sig_mask(:, col))
                sig_chans = find(sig_mask(:, col));
                topoplot(data_to_plot, chanlocs, 'maplimits', clim_Shared, 'colormap', jet, ...
                    'emarker2', {sig_chans, '.', 'k', 16, 2});
            else
                topoplot(data_to_plot, chanlocs, 'maplimits', clim_Shared, 'colormap', jet);
            end

            % Left row titles
            if col == 1
                text(-0.72, 0, row_titles{row}, 'FontSize', 16, 'FontWeight', 'bold', ...
                    'Rotation', 90, 'HorizontalAlignment', 'center');
            end

            % Row 3: Significant channel callouts (top) and time axis (bottom)
            if row == 3
                if any(sig_mask(:, col))
                    sig_ch_names = chan_labels(sig_mask(:, col));
                    n_sig = length(sig_ch_names);
                    if n_sig > 6
                        split_idx = ceil(n_sig / 2);
                        line1 = strjoin(sig_ch_names(1:split_idx), ', ');
                        line2 = strjoin(sig_ch_names(split_idx+1:end), ', ');
                        title({line1, line2}, 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
                    else
                        sig_str = strjoin(sig_ch_names, ', ');
                        title(sig_str, 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
                    end
                else
                    title('');
                end

                % Bottom time labels (using t_axis)
                text(0, -0.68, sprintf('%.2fs', t_axis(col)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
            else
                title('');
            end
        end
    end

    % --- 7. Colorbar Placement ---
    t.OuterPosition = [0, 0.04, 0.94, 0.94];
    cb_axes = [axes(fig), axes(fig), axes(fig)];
    set(cb_axes, 'Visible', 'off');

    colormap(cb_axes(1), jet); clim(cb_axes(1), clim_Shared);
    cb1 = colorbar(cb_axes(1), 'Position', [0.90 0.70 0.012 0.20]);
    cb1.FontSize = 14;
    ylabel(cb1, 'Normalized Amp (z)', 'FontSize', 14, 'FontWeight', 'bold');

    colormap(cb_axes(2), jet); clim(cb_axes(2), clim_Shared);
    cb2 = colorbar(cb_axes(2), 'Position', [0.90 0.38 0.012 0.20]);
    cb2.FontSize = 14;
    ylabel(cb2, 'Normalized Amp (z)', 'FontSize', 14, 'FontWeight', 'bold');

    colormap(cb_axes(3), jet); clim(cb_axes(3), clim_Shared);
    cb3 = colorbar(cb_axes(3), 'Position', [0.90 0.08 0.012 0.20]);
    cb3.FontSize = 14;
    ylabel(cb3, '\Delta z (Difference)', 'FontSize', 14, 'FontWeight', 'bold');

    sgtitle(sprintf('Group-Level Permutation Test: %s vs %s (%s)', condA_name, condB_name, state_name), ...
        'FontSize', 18, 'FontWeight', 'bold');

    % --- 8. Safe File Saving (Preventing Path Crash on '/') ---
    safe_state = strrep(strrep(strrep(strrep(strrep(state_name, ' ', '_'), '(', ''), ')', ''), '/', '_'), '\', '_');
    save_name = fullfile(output_dir, sprintf('Group_Power_Perm_%s_vs_%s_%s.png', ...
        strrep(condA_name, ' ', ''), strrep(condB_name, ' ', ''), safe_state));
    saveas(fig, save_name);
    close(fig);
end