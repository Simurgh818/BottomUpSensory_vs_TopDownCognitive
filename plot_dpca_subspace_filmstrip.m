function plot_dpca_subspace_filmstrip(avgA, avgB, W, k_opt, t_win, time_ms_eeg, window_size_ms, step_size_ms, bg_windows, cleanA, cleanB, all_channels_str, output_dir, subj_id)
    
    pc_labels = arrayfun(@(x) sprintf('dPC%d', x), 1:k_opt, 'UniformOutput', false);
    
    % --- 1. CALCULATE STATIC PRE-STIMULUS BASELINES (0.0 - 0.5s) ---
    idx_bg = find(time_ms_eeg >= bg_windows(1) - 1e-5 & time_ms_eeg <= bg_windows(2) + 1e-5);
    bg_pc_A = corrcoef((W * avgA(:, idx_bg))');
    bg_pc_B = corrcoef((W * avgB(:, idx_bg))');
    
    % --- 2. SLIDING WINDOW EVALUATION ---
    window_size_s = window_size_ms / 1000;
    step_size_s = step_size_ms / 1000;
    t_centers = (t_win(1) + step_size_s/2) : step_size_s : (t_win(2) - step_size_s/2 + 1e-5);
    num_windows = length(t_centers);
    
    sub_pc_A_all = nan(k_opt, k_opt, num_windows);
    sub_pc_B_all = nan(k_opt, k_opt, num_windows);
    
    % Storage for dPC traces and trajectory means
    trace_pc_A = cell(1, num_windows);
    trace_pc_B = cell(1, num_windows);
    time_chunks = cell(1, num_windows);
    
    dPC1_A = zeros(1, num_windows); dPC2_A = zeros(1, num_windows);
    dPC1_B = zeros(1, num_windows); dPC2_B = zeros(1, num_windows);
    
    for w = 1:num_windows
        wc = t_centers(w);
        idx_act = find(time_ms_eeg >= (wc - window_size_s/2 - 1e-5) & time_ms_eeg <= (wc + window_size_s/2 + 1e-5));
        
        act_A = avgA(:, idx_act);
        act_B = avgB(:, idx_act);
        time_chunks{w} = time_ms_eeg(idx_act);
        
        pc_A_act = W * act_A;
        pc_B_act = W * act_B;
        
        trace_pc_A{w} = pc_A_act;
        trace_pc_B{w} = pc_B_act;
        
        % Mean subspace states for trajectory plotting
        mean_pc_A = mean(pc_A_act, 2, 'omitnan');
        mean_pc_B = mean(pc_B_act, 2, 'omitnan');
        dPC1_A(w) = mean_pc_A(1); dPC2_A(w) = mean_pc_A(2);
        dPC1_B(w) = mean_pc_B(1); dPC2_B(w) = mean_pc_B(2);
        
        % Calculate and Subtract Static Baseline
        sub_pc_A_all(:,:,w) = corrcoef(pc_A_act') - bg_pc_A;
        sub_pc_B_all(:,:,w) = corrcoef(pc_B_act') - bg_pc_B;
    end
    
    % Concatenate all time chunks first, THEN extract the rows
    concat_A = cat(2, trace_pc_A{:});
    concat_B = cat(2, trace_pc_B{:});
    
    % Find global min/max for dPC1 and dPC2 traces to lock Y-axes
    all_dPC1 = [concat_A(1, :), concat_B(1, :)];
    all_dPC2 = [concat_A(2, :), concat_B(2, :)];
    pad1 = (max(all_dPC1) - min(all_dPC1)) * 0.1;
    pad2 = (max(all_dPC2) - min(all_dPC2)) * 0.1;
    ylim1 = [min(all_dPC1)-pad1, max(all_dPC1)+pad1]; if diff(ylim1)==0, ylim1=[-1 1]; end
    ylim2 = [min(all_dPC2)-pad2, max(all_dPC2)+pad2]; if diff(ylim2)==0, ylim2=[-1 1]; end
    
    %% --- 3. PLOTTING 5-ROW FIGURE ---
    c_lim = [-1 1]; 
    fig_width = max(1200, 220 * num_windows);
    fig = figure('Position', [50, 50, fig_width, 1100], 'Name', sprintf('%s vs %s dPCA Filmstrip', cleanA, cleanB), 'Visible', 'off');
    
    % Use a 5-row layout
    tiledlayout(5, num_windows, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    % ---------------------------------------------------------------------
    % ROW 1: Subtracted Condition A (Latent PCs)
    % ---------------------------------------------------------------------
    for w = 1:num_windows
        nexttile; imagesc(sub_pc_A_all(:,:,w)); colormap('jet'); clim(c_lim); axis square;
        if w == 1
            text(-0.75, k_opt/2, sprintf('%s\n(Shared dPCs)', cleanA), 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold');
            yticks(1:k_opt); yticklabels(pc_labels); set(gca, 'FontSize', 8);
        else
            yticks([]); 
        end
        xticks([]);
    end
    
    % ---------------------------------------------------------------------
    % ROW 2: Subtracted Condition B (Latent PCs) with TIME LABELS
    % ---------------------------------------------------------------------
    for w = 1:num_windows
        nexttile; imagesc(sub_pc_B_all(:,:,w)); colormap('jet'); clim(c_lim); axis square;
        if w == 1
            text(-0.75, k_opt/2, sprintf('%s\n(Shared dPCs)', cleanB), 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold');
            yticks(1:k_opt); yticklabels(pc_labels); set(gca, 'FontSize', 8);
        else
            yticks([]); 
        end
        xticks(1:k_opt); xticklabels(pc_labels); xtickangle(45); set(gca, 'FontSize', 8);
        
        % Time label centered below the X-ticks
        w_start = t_centers(w) - (step_size_s / 2);
        xlabel(sprintf('%.2fs', w_start), 'FontSize', 15, 'FontWeight', 'bold', 'Color', 'k');
    end
    
    % ---------------------------------------------------------------------
    % ROW 3: dPC1 TRACES (100 ms waveforms)
    % ---------------------------------------------------------------------
    for w = 1:num_windows
        nexttile; hold on;
        plot(time_chunks{w}, trace_pc_A{w}(1,:), 'Color', [0.2 0.6 0.8], 'LineWidth', 2);
        plot(time_chunks{w}, trace_pc_B{w}(1,:), 'Color', [0.8 0.4 0.2], 'LineWidth', 2);
        ylim(ylim1);
        
        if w == 1
            ylabel('dPC1 Amp.', 'FontSize', 12, 'FontWeight', 'bold');
        else
            yticks([]);
        end
        xticks([]); box on; set(gca, 'Color', [0.95 0.95 0.95]);
    end
    
    % ---------------------------------------------------------------------
    % ROW 4: dPC2 TRACES (100 ms waveforms)
    % ---------------------------------------------------------------------
    for w = 1:num_windows
        nexttile; hold on;
        plot(time_chunks{w}, trace_pc_A{w}(2,:), 'Color', [0.2 0.6 0.8], 'LineWidth', 2);
        plot(time_chunks{w}, trace_pc_B{w}(2,:), 'Color', [0.8 0.4 0.2], 'LineWidth', 2);
        ylim(ylim2);
        
        if w == 1
            ylabel('dPC2 Amp.', 'FontSize', 12, 'FontWeight', 'bold');
        else
            yticks([]);
        end
        xticks([]); box on; set(gca, 'Color', [0.95 0.95 0.95]);
    end
    
    % ---------------------------------------------------------------------
    % ROW 5: SPATIAL TRAJECTORY (dPC1 vs dPC2) - Spans entire bottom row
    % ---------------------------------------------------------------------
    nexttile([1, num_windows]); 
    hold on;
    
    % Plot the connecting lines and points
    plot(dPC1_A, dPC2_A, '-o', 'Color', [0.2 0.6 0.8], 'LineWidth', 2, 'MarkerFaceColor', [0.2 0.6 0.8], 'MarkerSize', 8, 'DisplayName', cleanA);
    plot(dPC1_B, dPC2_B, '-s', 'Color', [0.8 0.4 0.2], 'LineWidth', 2, 'MarkerFaceColor', [0.8 0.4 0.2], 'MarkerSize', 8, 'DisplayName', cleanB);
    
    % Annotate points with precise window start times
    for w = 1:num_windows
        w_start = t_centers(w) - (step_size_s / 2);
        text(dPC1_A(w), dPC2_A(w) + (pad2/4), sprintf(' %.2fs', w_start), 'Color', [0.2 0.6 0.8], 'FontSize', 10, 'FontWeight', 'bold');
        text(dPC1_B(w), dPC2_B(w) + (pad2/4), sprintf(' %.2fs', w_start), 'Color', [0.8 0.4 0.2], 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    xlabel('dPC 1 Amplitude', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('dPC 2 Amplitude', 'FontSize', 12, 'FontWeight', 'bold');
    title('Latent Subspace Trajectory (dPC1 vs dPC2)', 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 12);
    grid on; hold off;
    
    % Master formatting
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Label.String = 'Subtracted Correlation (\Delta r)';
    cb.Label.FontSize = 14;
    cb.Label.FontWeight = 'bold';
    
    sgtitle(sprintf('Latent Subspace Dynamics (0.0-0.5s Baseline Subtracted): %s vs %s', cleanA, cleanB), 'FontSize', 22, 'FontWeight', 'bold');
    
    save_file = fullfile(output_dir, sprintf('%s_%s_vs_%s_dPCA_Subspace_Filmstrip.png', subj_id, strrep(cleanA,' ',''), strrep(cleanB,' ','')));
    saveas(fig, save_file);
    close(fig);
end