function plot_within_subj_bands(mean_diff_matrix, time_axis, band_names, cleanA, cleanB, state_name, output_dir, subj_id)
    % mean_diff_matrix dimensions: [Bands x Time Windows]
    
    num_bands = length(band_names);
    fig = figure('Position', [150, 150, 1000, 800], 'Name', sprintf('%s Band Profile', state_name), 'Visible', 'off');
    tiledlayout(num_bands, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for b = 1:num_bands
        nexttile;
        plot(time_axis, mean_diff_matrix(b, :), 'k', 'LineWidth', 2);
        
        ylabel(sprintf('%s\n(|\\Delta r|)', upper(band_names{b})), 'FontWeight', 'bold', 'FontSize', 16);
        grid on;
        
        % Dynamic limits to handle 1/f drop-off
        % y_max = max(mean_diff_matrix(b, :)) * 1.2;
        % if y_max == 0 || isnan(y_max), y_max = 0.1; end
        ylim([0, 1]); % y_max
        xlim([time_axis(1), time_axis(end)]);
        
        if b == 1
            title(sprintf('Within-Subject Profile (%s): %s vs %s', state_name, cleanA, cleanB), 'FontSize', 16, 'FontWeight', 'bold');
        end
        if b == num_bands
            xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold');
        else
            xticklabels({}); % Hide x-axis labels on upper subplots
        end
    end
    
    % Save to the subject's main directory
    save_name = fullfile(output_dir, sprintf('%s_%s_vs_%s_BandProfile_%s.png', subj_id, strrep(cleanA,' ',''), strrep(cleanB,' ',''), strrep(state_name,' ','_')));
    saveas(fig, save_name);
    close(fig);
end