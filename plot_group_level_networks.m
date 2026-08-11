function plot_group_level_networks(grand_avg_net, p_values, window_centers, condA_name, condB_name, state_name, current_band, all_channels_str, chanlocs, output_dir)
    
    num_windows = length(window_centers);
    if num_windows == 0, return; end
    num_ch = size(grand_avg_net, 1);
    
    % --- FDR Correction Parameters ---
    q_fdr = 0.10; 
    sig_mask = zeros(size(p_values));
    
    % Apply Benjamini-Hochberg FDR correction per window
    for w = 1:num_windows
        p_win = p_values(:,:,w);
        
        % Extract only the upper triangle to avoid double-counting symmetric edges
        mask_tri = triu(true(size(p_win)), 1);
        p_tri = p_win(mask_tri);
        
        % Sort p-values and apply BH threshold
        [p_sort, ~] = sort(p_tri);
        m = length(p_sort);
        ranks = (1:m)';
        
        sig_idx = find(p_sort <= (ranks / m) * q_fdr, 1, 'last');
        
        if ~isempty(sig_idx)
            p_thresh = p_sort(sig_idx);
            % Create logical mask of edges surviving FDR
            sig_mask(:,:,w) = (p_win <= p_thresh) & mask_tri;
            % Make mask symmetric for plotting
            sig_mask(:,:,w) = sig_mask(:,:,w) | sig_mask(:,:,w)'; 
        end
    end
    
    % --- Extract Electrode Coordinates for Head Model ---
    if ~isempty(chanlocs)
        th = pi/180 * [chanlocs.theta];
        rd = [chanlocs.radius];
        x_node = rd .* sin(th);
        y_node = rd .* cos(th);
        rmax = max(0.5, max(rd)); 
    else
        % Fallback (Draws a circle on the circumference)
        th = linspace(0, 2*pi, num_ch+1); th(end)=[];
        x_node = 0.5 * cos(th);
        y_node = 0.5 * sin(th);
        rmax = 0.5;
    end
    
    % --- Visualization Setup ---
    c_lim_mat = [-1 1]; 
    cmap = jet(256);
    get_color = @(v) cmap(max(1, min(256, round((v + 1) / 2 * 255) + 1)), :);
    
    fig_width = max(1200, 250 * num_windows);
    fig1 = figure('Position', [100, 100, fig_width, 400], 'Name', sprintf('Group Level %s Band', current_band), 'Visible', 'off');
    
    tiledlayout(1, num_windows, 'TileSpacing', 'tight', 'Padding', 'tight');
    
    theta_circle = linspace(0, 2*pi, 100);
    draw_head = @() plot(rmax*cos(theta_circle), rmax*sin(theta_circle), 'k', 'LineWidth', 1.5);
    draw_nose = @() plot([rmax*0.1, 0, -rmax*0.1], [rmax, rmax*1.15, rmax], 'k', 'LineWidth', 1.5);
    
    % --- Plot FDR-Corrected Group Networks ---
    for w = 1:num_windows
        nexttile; hold on;
        
        draw_head(); draw_nose();
        scatter(x_node, y_node, 20, 'k', 'filled'); 
        
        % Extract the grand average network and apply the FDR mask
        mat_avg = grand_avg_net(:,:,w);
        mat_mask = sig_mask(:,:,w);
        
        [row, col] = find(mat_mask == 1);
        
        % Draw the surviving significant edges
        for i = 1:length(row)
            if row(i) > col(i) 
                v = mat_avg(row(i), col(i));
                % Line width scales with the absolute grand-average correlation difference
                plot([x_node(row(i)) x_node(col(i))], [y_node(row(i)) y_node(col(i))], ...
                    'Color', get_color(v), 'LineWidth', max(1, abs(v) * 5));
            end
        end
        
        axis equal; axis off;
        colormap(gca, 'jet'); clim(c_lim_mat); 
        
        title(sprintf('%.2fs', window_centers(w)), 'FontSize', 16, 'FontWeight', 'bold');
        
        if w == 1
            % Annotate the FDR significance threshold
            text(-rmax*1.5, 0, sprintf('FDR (q < %g)', q_fdr), 'Rotation', 90, ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 16);
        end
        if w == num_windows
            ax_master_color = gca;
        end
    end
    
    % --- Add Master Colorbar & Title ---
    cb = colorbar(ax_master_color);
    cb.Layout.Tile = 'east';
    cb.Label.String = 'Group Average \Delta r';
    cb.Label.FontSize = 16;
    cb.Label.FontWeight = 'bold';
    
    sgtitle(sprintf('[%s Band] Significant Group Networks: %s vs %s (%s)', ...
        upper(current_band), condA_name, condB_name, state_name), 'FontSize', 22, 'FontWeight', 'bold');
    
    % Save the Figure
    save_name = fullfile(output_dir, sprintf('GroupStat_%s_vs_%s_%s_%s.png', ...
        strrep(condA_name,' ',''), strrep(condB_name,' ',''), strrep(state_name,' ','_'), current_band));
    saveas(fig1, save_name);
    close(fig1);
end