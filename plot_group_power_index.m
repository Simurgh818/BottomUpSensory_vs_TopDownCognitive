function plot_group_power_index(pow_A_group, pow_B_group, t_axis, condA_name, condB_name, state_name, chanlocs, n_perms, alpha_level, output_dir)
    % PLOT_GROUP_POWER_INDEX Generates group-level Band Power Index (Alpha & Beta)
    % and Beta/Alpha Ratio Index topoplots with permutation testing.

    if nargin < 8 || isempty(n_perms), n_perms = 1000; end
    if nargin < 9 || isempty(alpha_level), alpha_level = 0.05; end
    if nargin < 10 || isempty(output_dir), output_dir = pwd; end
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    num_windows = length(t_axis);
    num_subjs = size(pow_A_group, 1);
    num_channels = size(pow_A_group{1, 1}, 1);
    num_tests = num_channels * num_windows;

    % Safe Channel Label Extraction
    if ~isempty(chanlocs) && isstruct(chanlocs) && isfield(chanlocs, 'labels')
        chan_labels = {chanlocs.labels};
    else
        chan_labels = arrayfun(@(x) sprintf('Ch%d', x), 1:num_channels, 'UniformOutput', false);
    end

    % Define Frequency Bands (Alpha = index 1, Beta = index 2)
    band_labels = {'ALPHA', 'BETA'};
    num_bands = 2;

    % =====================================================================
    % 1. COMPUTE SUBJECT TRIAL AVERAGES & INDICES
    % =====================================================================
    subj_powA = zeros(num_subjs, num_channels, num_windows, num_bands);
    subj_powB = zeros(num_subjs, num_channels, num_windows, num_bands);
    subj_idx  = zeros(num_subjs, num_channels, num_windows, num_bands);

    for s = 1:num_subjs
        for b = 1:num_bands
            datA = pow_A_group{s, b};
            datB = pow_B_group{s, b};

            % Trial average per subject: [Channels x Windows]
            if ndims(datA) == 3
                mA = mean(datA, 3, 'omitnan');
                mB = mean(datB, 3, 'omitnan');
            else
                mA = datA;
                mB = datB;
            end

            subj_powA(s, :, :, b) = mA;
            subj_powB(s, :, :, b) = mB;

            % Power Index = (Cond A - Cond B) / (Cond A + eps)
            subj_idx(s, :, :, b) = (mA - mB) ./ (mA + eps);
        end
    end

    % Beta/Alpha Ratio Index
    subj_ratioA = squeeze(subj_powA(:, :, :, 2)) ./ (squeeze(subj_powA(:, :, :, 1)) + eps);
    subj_ratioB = squeeze(subj_powB(:, :, :, 2)) ./ (squeeze(subj_powB(:, :, :, 1)) + eps);
    subj_ratio_idx = (subj_ratioA - subj_ratioB) ./ (subj_ratioA + eps);

    % Observed Grand Averages across subjects
    obs_band_idx = squeeze(mean(subj_idx, 1, 'omitnan')); % [Channels x Windows x Bands]
    obs_ratio_idx = squeeze(mean(subj_ratio_idx, 1, 'omitnan')); % [Channels x Windows]

    % =====================================================================
    % 2. MAX-STATISTIC PERMUTATION TESTING (Alpha, Beta, and Ratio)
    % =====================================================================
    sig_mask_bands = false(num_channels, num_windows, num_bands);
    sig_mask_ratio = false(num_channels, num_windows);

    % Run Permutations per Band
    for b = 1:num_bands
        max_null = zeros(n_perms, 1);
        for p = 1:n_perms
            pseudo_subjs = zeros(num_subjs, num_channels, num_windows);
            for s = 1:num_subjs
                mA = squeeze(subj_powA(s, :, :, b));
                mB = squeeze(subj_powB(s, :, :, b));
                if rand() > 0.5
                    pseudo_subjs(s, :, :) = (mA - mB) ./ (mA + eps);
                else
                    pseudo_subjs(s, :, :) = (mB - mA) ./ (mB + eps);
                end
            end
            pseudo_grp = squeeze(mean(pseudo_subjs, 1, 'omitnan'));
            max_null(p) = max(abs(pseudo_grp(:)));
        end
        thresh = prctile(max_null, (1 - alpha_level) * 100);
        sig_mask_bands(:, :, b) = abs(obs_band_idx(:, :, b)) > thresh;
    end

    % Run Permutations for Beta/Alpha Ratio
    max_null_ratio = zeros(n_perms, 1);
    for p = 1:n_perms
        pseudo_subjs_r = zeros(num_subjs, num_channels, num_windows);
        for s = 1:num_subjs
            rA = squeeze(subj_ratioA(s, :, :));
            rB = squeeze(subj_ratioB(s, :, :));
            if rand() > 0.5
                pseudo_subjs_r(s, :, :) = (rA - rB) ./ (rA + eps);
            else
                pseudo_subjs_r(s, :, :) = (rB - rA) ./ (rB + eps);
            end
        end
        pseudo_grp_r = squeeze(mean(pseudo_subjs_r, 1, 'omitnan'));
        max_null_ratio(p) = max(abs(pseudo_grp_r(:)));
    end
    thresh_ratio = prctile(max_null_ratio, (1 - alpha_level) * 100);
    sig_mask_ratio = abs(obs_ratio_idx) > thresh_ratio;

    % Dynamic Figure Dimensions (Prevents Squeeze on 21+ Windows)
    fig_width = max(1100, 180 * num_windows + 150);

    % =====================================================================
    % 3. FIGURE 1: BAND POWER INDEX (ALPHA & BETA - 2 ROWS)
    % =====================================================================
    fig1 = figure('Position', [50, 80, fig_width, 500], 'Name', 'Band Power Index', 'Visible', 'off');
    t1 = tiledlayout(num_bands, num_windows, 'TileSpacing', 'compact', 'Padding', 'compact');
    t1.OuterPosition = [0, 0.02, 0.88, 0.96]; % 88% width protects against overlap

    for b = 1:num_bands
        for col = 1:num_windows
            nexttile(t1);

            data_to_plot = obs_band_idx(:, col, b);
            curr_sig = sig_mask_bands(:, col, b);

            if any(curr_sig)
                sig_chans = find(curr_sig);
                topoplot(data_to_plot, chanlocs, 'maplimits', [-1 1], 'colormap', jet, ...
                    'emarker2', {sig_chans, '.', 'k', 15, 2});

                % Significant channel callout
                sig_ch_names = chan_labels(curr_sig);
                if length(sig_ch_names) > 5
                    split_idx = ceil(length(sig_ch_names)/2);
                    line1 = strjoin(sig_ch_names(1:split_idx), ', ');
                    line2 = strjoin(sig_ch_names(split_idx+1:end), ', ');
                    title({line1, line2}, 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
                else
                    title(strjoin(sig_ch_names, ', '), 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
                end
            else
                topoplot(data_to_plot, chanlocs, 'maplimits', [-1 1], 'colormap', jet);
                title('');
            end

            % Left-side Band Name
            if col == 1
                text(-0.75, 0, band_labels{b}, 'HorizontalAlignment', 'center', ...
                    'Rotation', 90, 'FontSize', 15, 'FontWeight', 'bold');
            end

            % Bottom-row time annotations
            if b == num_bands
                text(0, -0.68, sprintf('%.2fs', t_axis(col)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
            end
        end
    end

    % Dedicated Non-Overlapping Colorbar for Figure 1
    cb_ax1 = axes(fig1); set(cb_ax1, 'Visible', 'off');
    colormap(cb_ax1, jet); clim(cb_ax1, [-1 1]);
    cb1 = colorbar(cb_ax1, 'Position', [0.91, 0.25, 0.012, 0.50]);
    cb1.FontSize = 13;
    cb1.Ticks = [-1, -0.5, 0, 0.5, 1];
    ylabel(cb1, sprintf('Power Index\n(%s - %s)/%s', condA_name, condB_name, condA_name), ...
        'FontSize', 13, 'FontWeight', 'bold');

    title(t1, sprintf('Band Power Index: %s vs %s (%s)', condA_name, condB_name, state_name), ...
        'FontSize', 18, 'FontWeight', 'bold');

    % Save Figure 1
    safe_state = strrep(strrep(strrep(strrep(state_name, ' ', '_'), '(', ''), ')', ''), '/', '_');
    save1 = fullfile(output_dir, sprintf('Band_Power_Index_%s_vs_%s_%s.png', ...
        strrep(condA_name, ' ', ''), strrep(condB_name, ' ', ''), safe_state));
    saveas(fig1, save1);
    close(fig1);

    % =====================================================================
    % 4. FIGURE 2: BETA / ALPHA RATIO INDEX (1 ROW)
    % =====================================================================
    fig2 = figure('Position', [50, 150, fig_width, 360], 'Name', 'Beta/Alpha Ratio', 'Visible', 'off');
    t2 = tiledlayout(1, num_windows, 'TileSpacing', 'compact', 'Padding', 'compact');
    t2.OuterPosition = [0, 0.02, 0.88, 0.94];

    for col = 1:num_windows
        nexttile(t2);

        data_to_plot = obs_ratio_idx(:, col);
        curr_sig = sig_mask_ratio(:, col);

        if any(curr_sig)
            sig_chans = find(curr_sig);
            topoplot(data_to_plot, chanlocs, 'maplimits', [-1 1], 'colormap', jet, ...
                'emarker2', {sig_chans, '.', 'k', 15, 2});

            sig_ch_names = chan_labels(curr_sig);
            if length(sig_ch_names) > 5
                split_idx = ceil(length(sig_ch_names)/2);
                line1 = strjoin(sig_ch_names(1:split_idx), ', ');
                line2 = strjoin(sig_ch_names(split_idx+1:end), ', ');
                title({line1, line2}, 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
            else
                title(strjoin(sig_ch_names, ', '), 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
            end
        else
            topoplot(data_to_plot, chanlocs, 'maplimits', [-1 1], 'colormap', jet);
            title('');
        end

        % Time labels on bottom
        text(0, -0.68, sprintf('%.2fs', t_axis(col)), ...
            'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    end

    % Dedicated Non-Overlapping Colorbar for Figure 2
    cb_ax2 = axes(fig2); set(cb_ax2, 'Visible', 'off');
    colormap(cb_ax2, jet); clim(cb_ax2, [-1 1]);
    cb2 = colorbar(cb_ax2, 'Position', [0.91, 0.20, 0.012, 0.60]);
    cb2.FontSize = 13;
    cb2.Ticks = [-1, -0.5, 0, 0.5, 1];
    ylabel(cb2, sprintf('Ratio Index\n(%s - %s)/%s', condA_name, condB_name, condA_name), ...
        'FontSize', 13, 'FontWeight', 'bold');

    title(t2, sprintf('Beta/Alpha Ratio: %s vs %s (%s)', condA_name, condB_name, state_name), ...
        'FontSize', 18, 'FontWeight', 'bold');

    % Save Figure 2
    save2 = fullfile(output_dir, sprintf('BetaAlpha_Ratio_Index_%s_vs_%s_%s.png', ...
        strrep(condA_name, ' ', ''), strrep(condB_name, ' ', ''), safe_state));
    saveas(fig2, save2);
    close(fig2);
end