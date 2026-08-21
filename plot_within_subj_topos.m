function plot_within_subj_topos(sA_cell, sB_cell, time_axis, band_names, cleanA, cleanB, state_name, chanlocs, output_dir, subj_id)
    
    num_bands = length(band_names);
    num_windows = length(time_axis);
    
    fig = figure('Position', [50, 50, max(600, 270 * num_windows), 220 * num_bands], ...
        'Name', sprintf('%s Band Topos', state_name), 'Visible', 'off');
    tiledlayout(num_bands, num_windows, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for b = 1:num_bands
        for w = 1:num_windows
            nexttile;
            
            % Extract the Active minus Background matrices for Control (A) and Test (B)
            matA = sA_cell{b}(:,:,w);
            matB = sB_cell{b}(:,:,w);
            
            % 1. Calculate Nodal Degree (Hub Strength) ignoring the diagonal
            nodeA = sum(matA - diag(diag(matA)), 2, 'omitnan');
            nodeB = sum(matB - diag(diag(matB)), 2, 'omitnan');
            
            % 2. Calculate the Stabilized Index: (Control - Test) / Control
            % Added +0.05 to the absolute denominator to prevent division-by-zero artifacts
            node_index = (nodeA - nodeB) ./ (abs(nodeA) + 0.05);
            topoplot(node_index, chanlocs, 'numcontour', 0);
            
            % STRICT LIMITS: Locked to exactly -1 to 1
            clim([-1 1]); 
            colormap('jet');
            
            % Format annotations
            if w == 1
                % --- UPDATED: Moved band name to the LEFT side of the first column ---
                text(-0.75, 0, upper(band_names{b}), ...
                     'HorizontalAlignment', 'center', 'Rotation', 90, ...
                     'FontSize', 16, 'FontWeight', 'bold');
            end
            
            if b == num_bands
                % Add the time strictly to the bottom of the last row
                text(0, -0.65, sprintf('%.2fs', time_axis(w)), 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
            end
            
            if b == num_bands && w == num_windows
                cb = colorbar;
                cb.Label.String = sprintf('Nodal Index\n(%s - %s) / %s', cleanA, cleanB, cleanA);
                cb.Label.FontSize = 14;
                cb.Label.FontWeight = 'bold';
            end
        end
    end
    
    sgtitle(state_name, 'FontSize', 22, 'FontWeight', 'bold');
    
    save_name = fullfile(output_dir, sprintf('%s_%s_vs_%s_DeltaR_Index_Topos_%s.png', subj_id, strrep(cleanA,' ',''), strrep(cleanB,' ',''), strrep(state_name,' ','_')));
    saveas(fig, save_name);
    close(fig);
end