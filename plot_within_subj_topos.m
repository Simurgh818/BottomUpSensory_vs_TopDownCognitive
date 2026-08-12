function plot_within_subj_topos(diff_state_cell, time_axis, band_names, cleanA, cleanB, state_name, chanlocs, output_dir, subj_id)
    
    num_bands = length(band_names);
    num_windows = length(time_axis);

    fig = figure('Position', [50, 50, max(600, 270 * num_windows), 220 * num_bands], ...
        'Name', sprintf('%s Band Topos', state_name), 'Visible', 'off');
    tiledlayout(num_bands, num_windows, 'TileSpacing', 'compact', 'Padding', 'compact');

    for b = 1:num_bands
        diff_mat_all = diff_state_cell{b}; % The full sliding [32 x 32 x Windows] differential tensor

        for w = 1:num_windows
            nexttile;
            curr_mat = diff_mat_all(:,:,w);
            
            % Calculate Node Degree Difference: Sum of correlation differences per channel (ignoring self-correlation)
            node_diff = sum(curr_mat - diag(diag(curr_mat)), 2, 'omitnan');

            topoplot(node_diff, chanlocs, 'numcontour', 0);
            
            % STRICT LIMITS: Locked to exactly -1 to 1
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
                text(0, -0.65, sprintf('%.2fs', time_axis(w)), 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
            end
            
            if b == num_bands && w == num_windows
                cb = colorbar;
                cb.Label.String = sprintf('Nodal \\Delta r\n(+ %s  |  - %s)', cleanB, cleanA);
                cb.Label.FontSize = 14;
                cb.Label.FontWeight = 'bold';
            end
        end
    end
    
    sgtitle(sprintf('Within-Subject Nodal \\Delta r: %s vs %s (%s)', cleanA, cleanB, state_name), 'FontSize', 22, 'FontWeight', 'bold');
    
    save_name = fullfile(output_dir, sprintf('%s_%s_vs_%s_DeltaR_Topos_%s.png', subj_id, strrep(cleanA,' ',''), strrep(cleanB,' ',''), strrep(state_name,' ','_')));
    saveas(fig, save_name);
    close(fig);
end