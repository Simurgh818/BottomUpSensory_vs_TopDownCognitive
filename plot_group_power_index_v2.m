function plot_group_power_index_v2(pow_A_group, pow_B_group, t_axis, condA_name, condB_name, state_name, chanlocs, n_perms, alpha_level, output_dir)
    % PLOT_GROUP_POWER_INDEX_V2 Generates a 4-row group-level topoplot figure:
    %   Row 1: RAW Broadband Signal Power Index  (CondA - CondB) / CondA
    %   Row 2: ALPHA Band Power Index            (CondA - CondB) / CondA
    %   Row 3: BETA Band Power Index             (CondA - CondB) / CondA
    %   Row 4: BETA / ALPHA Ratio Index          (CondA - CondB) / CondA
    % Evaluated over [-0.2, 1.0 s] relative to tactile onset (t = 0 s).

    if nargin < 8 || isempty(n_perms), n_perms = 1000; end
    if nargin < 9 || isempty(alpha_level), alpha_level = 0.05; end
    if nargin < 10 || isempty(output_dir), output_dir = pwd; end
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    % =====================================================================
    % 1. ALIGN TIME AXIS TO TACTILE ONSET (t = 0 s) & FILTER TO [-0.2, 1.0 s]
    % =====================================================================
    t_axis = t_axis(:)';
    if length(t_axis) > 1
        step_s = mode(diff(t_axis));
    else
        step_s = 0.20;
    end
    half_win_s = step_s / 2;

    % If t_axis is in raw epoch seconds (0 to 3.5 s), tactile onset is 1.0 s.
    % If t_axis is already tactile-locked (time_s from -1.0 to 2.5 s), tactile onset is 0.0 s.
    if min(t_axis) >= 0 && max(t_axis) > 1.2
        tactile_onset = 1.0;
    else
        tactile_onset = 0.0;
    end

    % Convert window centers to tactile-relative start times
    t_rel_centers = t_axis - tactile_onset;
    t_rel_starts  = t_rel_centers - half_win_s;

    % Restrict windows strictly to [-0.2, 1.0 s] relative to tactile onset
    win_keep = (t_rel_centers >= -0.200 - 1e-4) & (t_rel_centers <= 1.000 + 1e-4);
    if any(win_keep) && ~all(win_keep)
        t_rel_centers = t_rel_centers(win_keep);
        t_rel_starts  = t_rel_starts(win_keep);
        for s = 1:size(pow_A_group, 1)
            for b = 1:size(pow_A_group, 2)
                if ndims(pow_A_group{s, b}) == 3
                    pow_A_group{s, b} = pow_A_group{s, b}(:, win_keep, :);
                    pow_B_group{s, b} = pow_B_group{s, b}(:, win_keep, :);
                else
                    pow_A_group{s, b} = pow_A_group{s, b}(:, win_keep);
                    pow_B_group{s, b} = pow_B_group{s, b}(:, win_keep);
                end
            end
        end
    end

    num_windows  = length(t_rel_centers);
    num_subjs    = size(pow_A_group, 1);
    num_channels = size(pow_A_group{1, 1}, 1);

    % Format bottom time labels in seconds relative to t = 0 tactile onset
    time_labels_sec = cell(1, num_windows);
    for col = 1:num_windows
        w_start_s = round(t_rel_starts(col) * 20) / 20;
        time_labels_sec{col} = sprintf('%+.2fs', w_start_s);
        time_labels_sec{col} = strrep(time_labels_sec{col}, '+0.00s', '0.00s');
        time_labels_sec{col} = strrep(time_labels_sec{col}, '-0.00s', '0.00s');
    end

    % Safe Channel Label Extraction
    if ~isempty(chanlocs) && isstruct(chanlocs) && isfield(chanlocs, 'labels')
        chan_labels = {chanlocs.labels};
    else
        chan_labels = arrayfun(@(x) sprintf('Ch%d', x), 1:num_channels, 'UniformOutput', false);
    end

    row_labels = {'RAW', 'ALPHA', 'BETA', 'BETA / ALPHA'};
    num_input_bands = size(pow_A_group, 2); % 3 bands: Col 1 = Raw, Col 2 = Alpha, Col 3 = Beta

    % =====================================================================
    % 2. COMPUTE SUBJECT TRIAL AVERAGES & NORMALIZED INDICES
    % =====================================================================
    subj_powA = zeros(num_subjs, num_channels, num_windows, num_input_bands);
    subj_powB = zeros(num_subjs, num_channels, num_windows, num_input_bands);
    subj_idx  = zeros(num_subjs, num_channels, num_windows, num_input_bands);

    for s = 1:num_subjs
        for b = 1:num_input_bands
            datA = pow_A_group{s, b};
            datB = pow_B_group{s, b};

            if ndims(datA) == 3
                mA = mean(datA, 3, 'omitnan');
                mB = mean(datB, 3, 'omitnan');
            else
                mA = datA;
                mB = datB;
            end

            subj_powA(s, :, :, b) = mA;
            subj_powB(s, :, :, b) = mB;

            % Normalized Index = (Cond A - Cond B) / (Cond A + eps)
            subj_idx(s, :, :, b) = (mA - mB) ./ (mA + eps);
        end
    end

    % Row 4: Beta/Alpha Ratio Index (Beta = Col 3, Alpha = Col 2)
    subj_ratioA    = squeeze(subj_powA(:, :, :, 3)) ./ (squeeze(subj_powA(:, :, :, 2)) + eps);
    subj_ratioB    = squeeze(subj_powB(:, :, :, 3)) ./ (squeeze(subj_powB(:, :, :, 2)) + eps);
    subj_ratio_idx = (subj_ratioA - subj_ratioB) ./ (subj_ratioA + eps);

    % Stack all 4 rows: [num_channels x num_windows x 4]
    obs_all_rows = zeros(num_channels, num_windows, 4);
    obs_all_rows(:, :, 1) = squeeze(mean(subj_idx(:, :, :, 1), 1, 'omitnan')); % Row 1: Raw
    obs_all_rows(:, :, 2) = squeeze(mean(subj_idx(:, :, :, 2), 1, 'omitnan')); % Row 2: Alpha
    obs_all_rows(:, :, 3) = squeeze(mean(subj_idx(:, :, :, 3), 1, 'omitnan')); % Row 3: Beta
    obs_all_rows(:, :, 4) = squeeze(mean(subj_ratio_idx, 1, 'omitnan'));       % Row 4: Beta/Alpha

    % =====================================================================
    % 3. MAX-STATISTIC PERMUTATION TESTING (ALL 4 ROWS)
    % =====================================================================
    sig_mask_all = false(num_channels, num_windows, 4);

    % Permutations for Raw (b=1), Alpha (b=2), and Beta (b=3)
    for b = 1:num_input_bands
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
            pseudo_grp  = squeeze(mean(pseudo_subjs, 1, 'omitnan'));
            max_null(p) = max(abs(pseudo_grp(:)));
        end
        thresh = prctile(max_null, (1 - alpha_level) * 100);
        sig_mask_all(:, :, b) = abs(obs_all_rows(:, :, b)) > thresh;
    end

    % Permutations for Beta/Alpha Ratio (Row 4)
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
        pseudo_grp_r      = squeeze(mean(pseudo_subjs_r, 1, 'omitnan'));
        max_null_ratio(p) = max(abs(pseudo_grp_r(:)));
    end
    thresh_ratio = prctile(max_null_ratio, (1 - alpha_level) * 100);
    sig_mask_all(:, :, 4) = abs(obs_all_rows(:, :, 4)) > thresh_ratio;

    % =====================================================================
    % 4. RENDER 4-ROW FIGURE (EXPANDED TOPOPLOTS, RIGHT-ALIGNED COLORBAR)
    % =====================================================================
    fig_width  = max(1600, 210 * num_windows + 340);
    fig_height = 1100;

    fig = figure('Position', [30, 30, fig_width, fig_height], ...
        'Name', sprintf('Figure 2: %s vs %s', condA_name, condB_name), ...
        'Color', 'w', 'Visible', 'off');

    % Canvas coordinate spacing:
    % Expanded usable_w to 0.74 so topoplots are larger and fill horizontal space
    left_margin = 0.05;
    usable_w    = 0.74; 
    cell_w      = usable_w / num_windows;
    plot_w      = cell_w * 0.90;
    plot_h      = 0.170;

    % Vertical positions for Rows 1, 2, 3, 4 with dedicated title headroom
    row_y_bots = [0.695, 0.490, 0.285, 0.055];

    for row = 1:4
        for col = 1:num_windows
            pos_x = left_margin + (col - 1) * cell_w;
            pos_y = row_y_bots(row);
            ax = axes('Position', [pos_x, pos_y, plot_w, plot_h]); %#ok<LAXES>
            axes(ax); %#ok<LAXES>

            data_to_plot = obs_all_rows(:, col, row);
            curr_sig     = sig_mask_all(:, col, row);

            if any(curr_sig)
                sig_chans = find(curr_sig);
                topoplot(data_to_plot, chanlocs, 'maplimits', [-1 1], 'colormap', jet, ...
                    'emarker2', {sig_chans, '.', 'k', 14, 2});

                sig_ch_names = chan_labels(curr_sig);
                if length(sig_ch_names) > 5
                    split_idx = ceil(length(sig_ch_names) / 2);
                    line1 = strjoin(sig_ch_names(1:split_idx), ', ');
                    line2 = strjoin(sig_ch_names(split_idx+1:end), ', ');
                    title({line1, line2}, 'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
                else
                    title(strjoin(sig_ch_names, ', '), 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
                end
            else
                topoplot(data_to_plot, chanlocs, 'maplimits', [-1 1], 'colormap', jet);
                title('');
            end

            % Row label on the far left of column 1
            if col == 1
                text(-0.78, 0, row_labels{row}, 'HorizontalAlignment', 'center', ...
                    'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold');
            end

            % Time label below each topoplot
            text(0, -0.68, time_labels_sec{col}, ...
                'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
        end
    end

    % ---------------------------------------------------------------------
    % 5. RIGHT-ALIGNED COLORBAR & NON-COLLIDING TEXT LABELS
    % ---------------------------------------------------------------------
    % Positioned near the right margin (x = 0.835) to eliminate dead white space
    cb_x = 0.835; 
    cb_y = 0.080; 
    cb_w = 0.016; 
    cb_h = 0.760;

    cb_ax = axes(fig, 'Position', [0, 0, 1, 1], 'Visible', 'off');
    colormap(cb_ax, jet); clim(cb_ax, [-1 1]);
    cb = colorbar(cb_ax, 'Position', [cb_x, cb_y, cb_w, cb_h]);
    cb.FontSize   = 13;
    cb.Ticks      = [-1, -0.5, 0, 0.5, 1];
    cb.TickLabels = {'-1.0', '-0.5', '0.0', '+0.5', '+1.0'};
    ylabel(cb, sprintf('Normalized Index: (%s - %s) / %s', condA_name, condB_name, condA_name), ...
        'FontSize', 14, 'FontWeight', 'bold');

    % Canvas annotation layer for clean text placement
    ax_ann = axes(fig, 'Position', [0, 0, 1, 1], 'Visible', 'off');
    hold(ax_ann, 'on');

    % Shifted +0.038 units to the right so labels NEVER overlap the "+1.0" or "-1.0" tick text
    label_x_pos = cb_x + cb_w + 0.038;

    text(ax_ann, label_x_pos, cb_y + cb_h, condA_name, ...
        'FontSize', 15, 'FontWeight', 'bold', 'Color', [0.85 0 0], 'VerticalAlignment', 'middle');

    text(ax_ann, label_x_pos, cb_y, condB_name, ...
        'FontSize', 15, 'FontWeight', 'bold', 'Color', [0 0.25 0.85], 'VerticalAlignment', 'middle');

    % Narrative legend block above the colorbar
    text(ax_ann, cb_x - 0.010, cb_y + cb_h + 0.065, sprintf('Red: Higher in %s', condA_name), ...
        'FontSize', 15, 'FontWeight', 'bold', 'Color', [0.85 0 0], 'VerticalAlignment', 'bottom');

    text(ax_ann, cb_x - 0.010, cb_y + cb_h + 0.035, sprintf('Blue: Higher in %s', condB_name), ...
        'FontSize', 15, 'FontWeight', 'bold', 'Color', [0 0.25 0.85], 'VerticalAlignment', 'bottom');

    % ---------------------------------------------------------------------
    % 6. CONCISE 2-ROW MAIN TITLE & SECTION HEADERS
    % ---------------------------------------------------------------------
    center_topos_x = left_margin + usable_w / 2;

    % Row 1 of Main Title (Primary Headline)
    text(ax_ann, center_topos_x, 0.960, ...
        sprintf('Sensor-Space Dynamics: %s vs. %s', condA_name, condB_name), ...
        'HorizontalAlignment', 'center', 'FontSize', 19, 'FontWeight', 'bold', 'Color', 'k');

    % Row 2 of Main Title (Window & Formula details)
    text(ax_ann, center_topos_x, 0.925, ...
        sprintf('Normalized Index: (%s - %s) / %s  [-0.20s to 1.00s, t = 0s Tactile Onset]', ...
        condA_name, condB_name, condA_name), ...
        'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);

    % Section 1 Sub-header (above Rows 1, 2, 3: Raw, Alpha, Beta)
    text(ax_ann, center_topos_x, 0.885, ...
        'Broadband & Band Power Index (Raw, Alpha, Beta)', ...
        'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);

    % Section 2 Sub-header (above Row 4: Beta/Alpha Ratio)
    text(ax_ann, center_topos_x, 0.245, ...
        'Beta / Alpha Ratio Index', ...
        'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);

    % ---------------------------------------------------------------------
    % 7. EXPORT 600 PPI HIGH-RES PNG & VECTOR SVG
    % ---------------------------------------------------------------------
    safe_A   = strrep(strrep(condA_name, ' ', ''), '.', '');
    safe_B   = strrep(strrep(condB_name, ' ', ''), '.', '');
    out_base = fullfile(output_dir, sprintf('Figure2_TopoIndex_4Rows_%s_vs_%s', safe_A, safe_B));
    
    try
        exportgraphics(fig, [out_base, '_600ppi.png'], 'Resolution', 600);
        saveas(fig, [out_base, '.svg']);
        fprintf('  -> Saved 600 ppi PNG: %s_600ppi.png\n', out_base);
        fprintf('  -> Saved Vector SVG: %s.svg\n', out_base);
    catch
        print(fig, [out_base, '_600ppi.png'], '-dpng', '-r600');
    end
    close(fig);
end