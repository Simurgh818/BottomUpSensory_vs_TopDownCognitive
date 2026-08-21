function plot_band_power_topos(trialsA_raw, trialsB_raw, t_win, time_ms_eeg, fs, win_sizes_ms, step_size_ms, condA, condB, active_state, bands, chanlocs, output_dir, subj_id)
    band_names = fieldnames(bands);
    num_bands = length(band_names);
    
    % Define exact time centers analytically to prevent drift
    step_size_s = step_size_ms / 1000;
    t_centers = t_win(1) : step_size_s : (t_win(2) + 1e-5);
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
    % 1. MAIN BAND POWER FIGURE
    % =====================================================================
    % Increased figure height (250 * num_bands + 100) to prevent label cutoff
    fig = figure('Position', [100, 100, max(600, 270 * num_windows), 250 * num_bands + 100], ...
        'Name', sprintf('%s Power Topos', active_state), 'Visible', 'off');
    tiledlayout(num_bands, num_windows, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    for b = 1:num_bands
        f_range = bands.(band_names{b});
        win_sz_sec = win_sizes_ms.(band_names{b}) / 1000; % Extract dynamic window size for FFT
        
        for w = 1:num_windows
            wc = t_centers(w);
            
            % Locate sample bounds, clipping to data edges if the window exceeds the epoch length
            idx_start = find(time_ms_eeg >= (wc - win_sz_sec/2), 1, 'first');
            if isempty(idx_start), idx_start = 1; end
            
            idx_end = find(time_ms_eeg <= (wc + win_sz_sec/2), 1, 'last');
            if isempty(idx_end), idx_end = length(time_ms_eeg); end
            
            % Isolate raw data chunk [Channels x Time x Trials]
            chunkA = trialsA_raw(:, idx_start:idx_end, :);
            chunkB = trialsB_raw(:, idx_start:idx_end, :);
            
            N = size(chunkA, 2);
            
            if N < 2
                powerA = zeros(size(chunkA, 1), 1);
                powerB = zeros(size(chunkB, 1), 1);
            else
                % 1. Apply Hanning window to prevent edge leakage prior to FFT
                h_win = hanning(N)'; % 1 x N
                chunkA_w = chunkA .* h_win; 
                chunkB_w = chunkB .* h_win; 
                
                % 2. Calculate FFT (Scaled by N to get proper amplitude)
                Y_A = fft(chunkA_w, [], 2) / N;
                Y_B = fft(chunkB_w, [], 2) / N;
                
                % 3. Extract Frequency Bins
                f_vec = (0:N-1) * (fs / N);
                f_idx = f_vec >= f_range(1) & f_vec <= f_range(2);
                
                % 4. Sum Power across targeted bins, then Average across Trials
                powA_trials = sum(abs(Y_A(:, f_idx, :)).^2, 2);
                powB_trials = sum(abs(Y_B(:, f_idx, :)).^2, 2);
                
                powerA = mean(powA_trials, 3, 'omitnan');
                powerB = mean(powB_trials, 3, 'omitnan');
            end
            
            % Save raw power to use later for ratios
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
                % --- UPDATED: Show ONLY the band name on the LEFT of the first column ---
                text(-0.75, 0, sprintf('%s\n(W=%dms)', upper(band_names{b}), win_sizes_ms.(band_names{b})), ...
                     'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold');
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
    
    % =====================================================================
    % 2. NEW BETA / ALPHA RATIO FIGURE
    % =====================================================================
    alpha_idx = find(strcmpi(band_names, 'alpha'));
    beta_idx  = find(strcmpi(band_names, 'beta'));
    
    if ~isempty(alpha_idx) && ~isempty(beta_idx)
        % Increased figure height to 450 to prevent label cutoff
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
            
            % Compute the Beta / Alpha Ratio
            ratioA = powA_beta ./ (powA_alpha + eps);
            ratioB = powB_beta ./ (powB_alpha + eps);
            
            % Calculate Ratio Contrast Index = (Cond A - Cond B) / Cond A
            ratio_index = (ratioA - ratioB) ./ (ratioA + eps);
            
            nexttile;
            topoplot(ratio_index, chanlocs, 'numcontour', 0);
            clim([-1 1]); 
            colormap('jet');
            
            % --- UPDATED: Removed the redundant 'BETA/ALPHA' title block completely ---
            
            text(0, -0.65, sprintf('%.2fs', wc), 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
            
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