function plot_band_power_topos(trialsA_raw, trialsB_raw, t_win, time_ms_eeg, fs, window_size_ms, step_size_ms, condA, condB, active_state, bands, chanlocs, output_dir, subj_id)

    band_names = fieldnames(bands);
    num_bands = length(band_names);

    % Configure sliding windows precisely aligned with the network connectivity
    win_samples = round((window_size_ms / 1000) * fs);
    step_samples = round((step_size_ms / 1000) * fs);

    idx_start = find(time_ms_eeg >= t_win(1), 1, 'first');
    idx_end = find(time_ms_eeg <= t_win(2), 1, 'last');
    time_chunk = time_ms_eeg(idx_start:idx_end);

    start_idx = 1 : step_samples : (length(time_chunk) - win_samples + 1);
    num_windows = length(start_idx);

    get_clean = @(c) strrep(strrep(strrep(strrep(strrep(strrep(c, ...
        'BLA', 'Auditory'), 'BLT', 'Tactile'), 'P1', 'Cued'), ...
        'P2', 'Unpred.'), 'P3', 'Rand. Cued'), '_', ' ');
    cleanA = get_clean(condA);
    cleanB = get_clean(condB);

    fig = figure('Position', [100, 100, max(600, 270 * num_windows), 220 * num_bands], ...
        'Name', sprintf('%s Power Topos', active_state), 'Visible', 'off');
    tiledlayout(num_bands, num_windows, 'TileSpacing', 'compact', 'Padding', 'compact');

    for b = 1:num_bands
        f_range = bands.(band_names{b});
        
        filtA = filter_trials_band(trialsA_raw, f_range, fs);
        filtB = filter_trials_band(trialsB_raw, f_range, fs);

        for w = 1:num_windows
            w_s = idx_start + start_idx(w) - 1;
            w_e = w_s + win_samples - 1;
            wc = mean(time_ms_eeg(w_s:w_e));

            % Calculate Time-Domain Power (Variance) and average across trials
            powerA = mean(squeeze(var(filtA(:, w_s:w_e, :), 0, 2, 'omitnan')), 2, 'omitnan');
            powerB = mean(squeeze(var(filtB(:, w_s:w_e, :), 0, 2, 'omitnan')), 2, 'omitnan');

            % Calculate Power Index = (Cond A - Cond B) / Cond A
            power_index = (powerA - powerB) ./ (powerA + eps);

            nexttile;
            topoplot(power_index, chanlocs, 'numcontour', 0);
            clim([-1 1]); 
            colormap('jet');

            % Format annotations
            if w == 1
                % Show ONLY the band name on the first column
                t = title(sprintf('%s', upper(band_names{b})), 'FontSize', 16, 'FontWeight', 'bold');
                t.Position(1) = t.Position(1) - 0.4; 
            end
            
            if b == num_bands
                % Add the time strictly to the bottom of the last row
                text(0, -0.65, sprintf('%.2fs', wc), 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
            end
            
            if b == num_bands && w == num_windows
                cb = colorbar;
                cb.Label.String = sprintf('Power Index\n(%s - %s)/%s', cleanA, cleanB, cleanA);
                cb.Label.FontSize = 14;
                cb.Label.FontWeight = 'bold';
            end
        end
    end
    
    sgtitle(sprintf('Band Power Index: %s vs %s (%s)', cleanA, cleanB, active_state), 'FontSize', 22, 'FontWeight', 'bold');
    
    save_name = fullfile(output_dir, sprintf('%s_%s_vs_%s_PowerTopos_%s.png', subj_id, condA, condB, strrep(active_state,' ','_')));
    saveas(fig, save_name);
    close(fig);
end