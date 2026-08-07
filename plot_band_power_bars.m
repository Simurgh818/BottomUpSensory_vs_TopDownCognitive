function plot_band_power_bars(trialsA_raw, trialsB_raw, t_win, time_ms_eeg, fs, condA, condB, active_state, bands, all_channels_str, output_dir, subj_id)
    
    % Get exact time indices for the window
    idx_start = find(time_ms_eeg >= t_win(1), 1, 'first');
    idx_end = find(time_ms_eeg <= t_win(2), 1, 'last');

    band_names = fieldnames(bands);
    num_bands = length(band_names);
    num_ch = size(trialsA_raw, 1);

    % Clean condition names for the legend
    get_clean = @(c) strrep(strrep(strrep(strrep(strrep(strrep(c, ...
        'BLA', 'Auditory'), 'BLT', 'Tactile'), 'P1', 'Cued'), ...
        'P2', 'Unpred.'), 'P3', 'Rand. Cued'), '_', ' ');
    cleanA = get_clean(condA);
    cleanB = get_clean(condB);

    % Create the 5-Row Figure
    fig = figure('Position', [100, 100, 1600, 1000], 'Name', sprintf('%s Band Power', active_state), 'Visible', 'off');
    tiledlayout(num_bands, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    for b = 1:num_bands
        f_range = bands.(band_names{b});

        % Filter the raw trials into the current frequency band
        filtA = filter_trials_band(trialsA_raw, f_range, fs);
        filtB = filter_trials_band(trialsB_raw, f_range, fs);

        % Calculate Power (Variance over the time dimension) for each trial
        % Resulting size: [Channels x Trials]
        powerA = squeeze(var(filtA(:, idx_start:idx_end, :), 0, 2, 'omitnan'));
        powerB = squeeze(var(filtB(:, idx_start:idx_end, :), 0, 2, 'omitnan'));

        % Calculate Mean and Standard Error (SEM) across trials
        meanA = mean(powerA, 2, 'omitnan');
        semA = std(powerA, 0, 2, 'omitnan') ./ sqrt(size(powerA, 2));

        meanB = mean(powerB, 2, 'omitnan');
        semB = std(powerB, 0, 2, 'omitnan') ./ sqrt(size(powerB, 2));

        % Plotting the grouped bar chart
        nexttile; hold on;
        b_bar = bar([meanA, meanB], 'grouped');
        b_bar(1).FaceColor = [0.2 0.6 0.8]; % Blueish
        b_bar(2).FaceColor = [0.8 0.4 0.2]; % Orangeish

        % Add SEM Error bars aligned to the center of each grouped bar
        x_posA = b_bar(1).XEndPoints;
        x_posB = b_bar(2).XEndPoints;
        errorbar(x_posA, meanA, semA, 'k', 'linestyle', 'none', 'CapSize', 2, 'LineWidth', 1);
        errorbar(x_posB, meanB, semB, 'k', 'linestyle', 'none', 'CapSize', 2, 'LineWidth', 1);

        ylabel(sprintf('%s\nPower (\\muV^2)', upper(band_names{b})), 'FontWeight', 'bold', 'FontSize', 12);
        xticks(1:num_ch);
        xlim([0.5, num_ch+0.5]);
        grid on;

        if b == num_bands
            xticklabels(all_channels_str);
            xtickangle(90);
            xlabel('Channels', 'FontWeight', 'bold');
        else
            xticklabels({});
        end

        if b == 1
            title(sprintf('Band Power (%s): %s vs %s', active_state, cleanA, cleanB), 'FontSize', 18, 'FontWeight', 'bold');
            legend({cleanA, cleanB}, 'Location', 'northeast', 'FontSize', 12);
        end
    end
    
    % Save the Figure
    save_name = fullfile(output_dir, sprintf('%s_%s_vs_%s_PowerBars_%s.png', subj_id, condA, condB, strrep(active_state,' ','_')));
    saveas(fig, save_name);
    close(fig);
end