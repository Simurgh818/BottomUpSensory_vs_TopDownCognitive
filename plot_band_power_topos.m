function plot_band_power_topos(trialsA_raw, trialsB_raw, t_win, time_ms_eeg, fs, ~, step_size_ms, condA, condB, active_state, bands, chanlocs, output_dir, subj_id)
    band_names = fieldnames(bands);
    num_bands = length(band_names);
    
    % Define exact time centers analytically to guarantee perfect 100ms contiguous blocks
    step_size_s = step_size_ms / 1000;
    t_centers = (t_win(1) + step_size_s/2) : step_size_s : (t_win(2) - step_size_s/2 + 1e-5);
    num_windows = length(t_centers);
    
    get_clean = @(c) strrep(strrep(strrep(strrep(strrep(strrep(c, ...
        'BLA', 'Auditory'), 'BLT', 'Tactile'), 'P1', 'Cued'), ...
        'P2', 'Unpred.'), 'P3', 'Rand. Cued'), '_', ' ');
    cleanA = get_clean(condA);
    cleanB = get_clean(condB);
    
    % Preallocate storage for the Beta/Alpha ratio calculation
    store_powA = cell(num_bands, num_windows);
    store_powB = cell(num_bands, num_windows);
    
    % =====================================================================
    % 0. WHOLE-TRIAL FILTERING & HILBERT POWER ENVELOPE (NEW METHOD)
    % =====================================================================
    fprintf('      -> Extracting whole-trial instantaneous power envelopes...\n');
    full_powA = zeros(size(trialsA_raw, 1), size(trialsA_raw, 2), num_bands);
    full_powB = zeros(size(trialsB_raw, 1), size(trialsB_raw, 2), num_bands);
    
    for b = 1:num_bands
        f_range = bands.(band_names{b});
        
        % 1. Filter the entire 3.5s epoch
        filtA = filter_trials_band(trialsA_raw, f_range, fs);
        filtB = filter_trials_band(trialsB_raw, f_range, fs);
        
        % 2. Extract instantaneous power via Hilbert envelope & square it
        powA_trials = zeros(size(filtA));
        powB_trials = zeros(size(filtB));
        
        % Hilbert must be computed column-wise (per trial, per channel)
        for tr = 1:size(filtA, 3)
            powA_trials(:,:,tr) = abs(hilbert(filtA(:,:,tr)')').^2;
        end
        for tr = 1:size(filtB, 3)
            powB_trials(:,:,tr) = abs(hilbert(filtB(:,:,tr)')').^2;
        end
        
        % 3. Average across trials to get stable continuous power time-series
        full_powA(:,:,b) = mean(powA_trials, 3, 'omitnan');
        full_powB(:,:,b) = mean(powB_trials, 3, 'omitnan');
    end

    % =====================================================================
    % 1. MAIN BAND POWER FIGURE
    % =====================================================================
    fig = figure('Position', [100, 100, max(600, 270 * num_windows), 250 * num_bands + 100], ...
        'Name', sprintf('%s Power Topos', active_state), 'Visible', 'off');
    tiledlayout(num_bands, num_windows, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    for b = 1:num_bands
        for w = 1:num_windows
            wc = t_centers(w);
            
            % Exact bounds for the 100ms visual block
            idx_start = find(time_ms_eeg >= (wc - step_size_s/2 - 1e-5), 1, 'first');
            idx_end   = find(time_ms_eeg <= (wc + step_size_s/2 + 1e-5), 1, 'last');
            
            % Average the high-resolution instantaneous power inside this 100ms window
            powerA = mean(full_powA(:, idx_start:idx_end, b), 2, 'omitnan');
            powerB = mean(full_powB(:, idx_start:idx_end, b), 2, 'omitnan');
            
            % Save power to use later for ratios
            store_powA{b, w} = powerA;
            store_powB{b, w} = powerB;
            
            % Calculate Power Index = (Cond A - Cond B) / Cond A
            power_index = (powerA - powerB) ./ (powerA + eps);
            
            nexttile;
            topoplot(power_index, chanlocs, 'numcontour', 0);
            clim([-1 1]); 
            colormap('jet');
            
            % Format annotations
            if w == 1
                text(-0.75, 0, sprintf('%s', upper(band_names{b})), ...
                     'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold');
            end
            
            if b == num_bands
                % Calculate and print the START of the window (-100ms, 0ms, etc.)
                w_start = wc - (step_size_s / 2);
                text(0, -0.65, sprintf('%.2fs', w_start), 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
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
    
    % =====================================================================
    % 2. NEW BETA / ALPHA RATIO FIGURE
    % =====================================================================
    alpha_idx = find(strcmpi(band_names, 'alpha'));
    beta_idx  = find(strcmpi(band_names, 'beta'));
    
    if ~isempty(alpha_idx) && ~isempty(beta_idx)
        fig_ratio = figure('Position', [150, 150, max(600, 270 * num_windows), 450], ...
            'Name', sprintf('%s Beta/Alpha Ratio', active_state), 'Visible', 'off');
        tiledlayout(1, num_windows, 'TileSpacing', 'compact', 'Padding', 'normal');
        
        for w = 1:num_windows
            wc = t_centers(w);
            
            % Retrieve the raw band powers computed in the previous loop
            powA_alpha = store_powA{alpha_idx, w};
            powA_beta  = store_powA{beta_idx, w};
            powB_alpha = store_powB{alpha_idx, w};
            powB_beta  = store_powB{beta_idx, w};
            
            ratioA = powA_beta ./ (powA_alpha + eps);
            ratioB = powB_beta ./ (powB_alpha + eps);
            
            ratio_index = (ratioA - ratioB) ./ (ratioA + eps);
            
            nexttile;
            topoplot(ratio_index, chanlocs, 'numcontour', 0);
            clim([-1 1]); 
            colormap('jet');
            
            w_start = wc - (step_size_s / 2);
            text(0, -0.65, sprintf('%.2fs', w_start), 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
            
            if w == num_windows
                cb = colorbar;
                cb.Label.String = sprintf('Ratio Index\n(%s - %s)/%s', cleanA, cleanB, cleanA);
                cb.Label.FontSize = 14;
                cb.Label.FontWeight = 'bold';
            end
        end
        
        sgtitle(sprintf('Beta/Alpha Ratio: %s vs %s (%s)', cleanA, cleanB, active_state), 'FontSize', 22, 'FontWeight', 'bold');
        save_name_ratio = fullfile(output_dir, sprintf('%s_%s_vs_%s_BetaAlphaRatio_%s.png', subj_id, condA, condB, strrep(active_state,' ','_')));
        saveas(fig_ratio, save_name_ratio);
        close(fig_ratio);
    end
end