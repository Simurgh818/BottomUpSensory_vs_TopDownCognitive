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

    if min(t_axis) >= 0 && max(t_axis) > 1.2
        tactile_onset = 1.0;
    else
        tactile_onset = 0.0;
    end

    t_rel_centers = t_axis - tactile_onset;
    t_rel_starts  = t_rel_centers - half_win_s;

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

    time_labels_sec = cell(1, num_windows);
    for col = 1:num_windows
        w_start_s = round(t_rel_starts(col) * 20) / 20;
        time_labels_sec{col} = sprintf('%+.2fs', w_start_s);
        time_labels_sec{col} = strrep(time_labels_sec{col}, '+0.00s', '0.00s');
        time_labels_sec{col} = strrep(time_labels_sec{col}, '-0.00s', '0.00s');
    end

    if ~isempty(chanlocs) && isstruct(chanlocs) && isfield(chanlocs, 'labels')
        chan_labels = {chanlocs.labels};
    else
        chan_labels = arrayfun(@(x) sprintf('Ch%d', x), 1:num_channels, 'UniformOutput', false);
    end

    row_labels = {'RAW', 'ALPHA', 'BETA', 'BETA / ALPHA'};
    num_input_bands = size(pow_A_group, 2);

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
            subj_idx(s, :, :, b) = (mA - mB) ./ (mA + eps);
        end
    end

    subj_ratioA    = squeeze(subj_powA(:, :, :, 3)) ./ (squeeze(subj_powA(:, :, :, 2)) + eps);
    subj_ratioB    = squeeze(subj_powB(:, :, :, 3)) ./ (squeeze(subj_powB(:, :, :, 2)) + eps);
    subj_ratio_idx = (subj_ratioA - subj_ratioB) ./ (subj_ratioA + eps);

    obs_all_rows = zeros(num_channels, num_windows, 4);
    obs_all_rows(:, :, 1) = squeeze(mean(subj_idx(:, :, :, 1), 1, 'omitnan'));
    obs_all_rows(:, :, 2) = squeeze(mean(subj_idx(:, :, :, 2), 1, 'omitnan'));
    obs_all_rows(:, :, 3) = squeeze(mean(subj_idx(:, :, :, 3), 1, 'omitnan'));
    obs_all_rows(:, :, 4) = squeeze(mean(subj_ratio_idx, 1, 'omitnan'));

    % =====================================================================
    % 3. STATISTICAL TESTING ACROSS CHANNELS & WINDOWS (WILCOXON + FDR)
    % =====================================================================
    sig_mask_all = false(num_channels, num_windows, 4);
    
    for row = 1:4
        if row <= 3
            subj_mat = squeeze(subj_idx(:, :, :, row));
        else
            subj_mat = subj_ratio_idx;
        end

        p_vals = zeros(num_channels, num_windows);

        for ch = 1:num_channels
            for w = 1:num_windows
                dist = squeeze(subj_mat(:, ch, w));
                if sum(~isnan(dist)) > 1 && any(dist ~= 0)
                    p_vals(ch, w) = signrank(dist);
                else
                    p_vals(ch, w) = 1.0;
                end
            end
        end
        
        for w = 1:num_windows
            p_win = p_vals(:, w);
            [p_sorted, ~] = sort(p_win, 'ascend');
            fdr_thresh_vec = (1:num_channels)' * (alpha_level / num_channels);
            max_idx = find(p_sorted <= fdr_thresh_vec, 1, 'last');
            if ~isempty(max_idx)
                p_crit = p_sorted(max_idx);
                sig_mask_all(:, w, row) = (p_win <= p_crit);
            end
        end
    end

    % =====================================================================
    % 4. RENDER 4-ROW FIGURE (EXPANDED HEADROOM & ADAPTIVE CHANNEL LABELS)
    % =====================================================================
    fig_width  = max(1750, 220 * num_windows + 430);
    fig_height = 1100;

    fig = figure('Position', [30, 30, fig_width, fig_height], ...
        'Name', sprintf('Figure 2: %s vs %s', condA_name, condB_name), ...
        'Color', 'w', 'Visible', 'off');

    left_margin = 0.045;
    usable_w    = 0.730; 
    cell_w      = usable_w / num_windows;
    plot_w      = cell_w * 0.90;
    plot_h      = 0.170;

    % Adjusted vertical spacing: increased gap between Row 2 (ALPHA) and Row 3 (BETA) from 38.5px to 55px
    row_y_bots = [0.705, 0.505, 0.285, 0.055];

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
                    'emarker2', {sig_chans, '.', 'k', 28, 2});

                sig_ch_names = chan_labels(curr_sig);
                n_sig = length(sig_ch_names);
                
                % --- ADAPTIVE 3-LINE MAXIMUM WITH WIDE HORIZONTAL PACKING ---
                max_lines = 3; % Strictly limit to at most 3 vertical lines
                
                if n_sig <= 4
                    % 1 Line: Large font
                    lines_cell = {strjoin(sig_ch_names, ', ')};
                    f_size = 11.0;
                elseif n_sig <= 8
                    % 2 Lines: Balanced font (~3-4 channels per line)
                    ch_per_line = ceil(n_sig / 2);
                    lines_cell = {strjoin(sig_ch_names(1:ch_per_line), ', '); ...
                                  strjoin(sig_ch_names(ch_per_line+1:end), ', ')};
                    f_size = 9.5;
                else
                    % 3 Lines: Increased width (4-5 channels per line) with compact font
                    ch_per_line = ceil(n_sig / max_lines);
                    n_lines_act = ceil(n_sig / ch_per_line);
                    lines_cell = cell(n_lines_act, 1);
                    for l_i = 1:n_lines_act
                        i_st = (l_i - 1) * ch_per_line + 1;
                        i_en = min(l_i * ch_per_line, n_sig);
                        lines_cell{l_i} = strjoin(sig_ch_names(i_st:i_en), ', ');
                    end
                    f_size = 8.0;
                end
                
                % Condition: Scale font size smaller if any individual line exceeds 24 characters
                max_chars = max(cellfun(@length, lines_cell));
                if max_chars > 24
                    f_size = max(7.0, f_size - 1.0);
                end
                
                title(lines_cell, 'FontSize', f_size, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
            else
                topoplot(data_to_plot, chanlocs, 'maplimits', [-1 1], 'colormap', jet);
                title('');
            end

            % Row label on the far left of column 1
            if col == 1
                text(-0.78, 0, row_labels{row}, 'HorizontalAlignment', 'center', ...
                    'Rotation', 90, 'FontSize', 22, 'FontWeight', 'bold');
            end

            % Time label strictly on the bottom row (Row 4)
            if row == 4
                text(0, -0.68, time_labels_sec{col}, ...
                    'HorizontalAlignment', 'center', 'FontSize', 21, 'FontWeight', 'bold');
            end
        end
    end

    % ---------------------------------------------------------------------
    % 5. RIGHT-ALIGNED COLORBAR & NON-COLLIDING TEXT LABELS
    % ---------------------------------------------------------------------
    cb_x = 0.818; 
    cb_y = 0.080; 
    cb_w = 0.016; 
    cb_h = 0.760;

    cb_ax = axes(fig, 'Position', [0, 0, 1, 1], 'Visible', 'off');
    colormap(cb_ax, jet); clim(cb_ax, [-1 1]);
    cb = colorbar(cb_ax, 'Position', [cb_x, cb_y, cb_w, cb_h]);
    cb.FontSize   = 19;
    cb.Ticks      = [-1, -0.5, 0, 0.5, 1];
    cb.TickLabels = {'-1.0', '-0.5', '0.0', '+0.5', '+1.0'};
    ylabel(cb, sprintf('Normalized Index: (%s - %s) / %s', condA_name, condB_name, condA_name), ...
        'FontSize', 22, 'FontWeight', 'bold');

    ax_ann = axes(fig, 'Position', [0, 0, 1, 1], 'Visible', 'off');
    hold(ax_ann, 'on');

    % Top label (+1.0 Cued): shifted right to completely clear '+1.0'
    label_x_top = cb_x + cb_w + 0.052;
    text(ax_ann, label_x_top, cb_y + cb_h, condA_name, ...
        'FontSize', 21, 'FontWeight', 'bold', 'Color', [0.85 0 0], 'VerticalAlignment', 'middle');

    % Bottom label (-1.0 Unpred 2000): placed at +0.038 so it never clips on the right edge
    label_x_bot = cb_x + cb_w + 0.038;
    text(ax_ann, label_x_bot, cb_y, condB_name, ...
        'FontSize', 21, 'FontWeight', 'bold', 'Color', [0 0.25 0.85], 'VerticalAlignment', 'middle');

    % Top Legend: shifted leftward to x = 0.748 and sized at 16 pt to prevent clipping
    leg_x = cb_x - 0.050;
    text(ax_ann, leg_x, cb_y + cb_h + 0.065, sprintf('Red: Higher in %s', condA_name), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.85 0 0], 'VerticalAlignment', 'bottom');

    text(ax_ann, leg_x, cb_y + cb_h + 0.035, sprintf('Blue: Higher in %s', condB_name), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', [0 0.25 0.85], 'VerticalAlignment', 'bottom');

    % ---------------------------------------------------------------------
    % 6. MAIN TITLE & SECTION HEADERS
    % ---------------------------------------------------------------------
    center_topos_x = left_margin + usable_w / 2;

    text(ax_ann, center_topos_x, 0.965, ...
        sprintf('Sensor-Space Dynamics: %s vs. %s', condA_name, condB_name), ...
        'HorizontalAlignment', 'center', 'FontSize', 27, 'FontWeight', 'bold', 'Color', 'k');

    text(ax_ann, center_topos_x, 0.930, ...
        sprintf('Normalized Index: (%s - %s) / %s  [-0.20s to 1.00s, t = 0s Tactile Onset]', ...
        condA_name, condB_name, condA_name), ...
        'HorizontalAlignment', 'center', 'FontSize', 21, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);

    text(ax_ann, center_topos_x, 0.890, ...
        'Broadband & Band Power Index (Raw, Alpha, Beta)', ...
        'HorizontalAlignment', 'center', 'FontSize', 22, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);

    text(ax_ann, center_topos_x, 0.245, ...
        'Beta / Alpha Ratio Index', ...
        'HorizontalAlignment', 'center', 'FontSize', 22, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);

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
end