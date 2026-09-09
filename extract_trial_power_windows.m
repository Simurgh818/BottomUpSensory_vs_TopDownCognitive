function [powA_all, powB_all, t_centers] = extract_trial_power_windows(trialsA_raw, trialsB_raw, t_win, time_ms_eeg, fs, step_size_ms, bands)
    
    band_names = fieldnames(bands);
    num_bands = length(band_names);
    
    step_size_s = step_size_ms / 1000;
    t_centers = (t_win(1) + step_size_s/2) : step_size_s : (t_win(2) - step_size_s/2 + 1e-5);
    num_windows = length(t_centers);
    
    nA = size(trialsA_raw, 3);
    nB = size(trialsB_raw, 3);
    num_channels = size(trialsA_raw, 1);
    
    powA_all = zeros(num_channels, num_windows, nA, num_bands);
    powB_all = zeros(num_channels, num_windows, nB, num_bands);
    
    for b = 1:num_bands
        f_range = bands.(band_names{b});
        
        filtA = filter_trials_band(trialsA_raw, f_range, fs);
        filtB = filter_trials_band(trialsB_raw, f_range, fs);
        
        powA_inst = zeros(size(filtA));
        powB_inst = zeros(size(filtB));
        
        for tr = 1:nA, powA_inst(:,:,tr) = abs(hilbert(filtA(:,:,tr)')').^2; end
        for tr = 1:nB, powB_inst(:,:,tr) = abs(hilbert(filtB(:,:,tr)')').^2; end
        
        for w = 1:num_windows
            wc = t_centers(w);
            idx_start = find(time_ms_eeg >= (wc - step_size_s/2 - 1e-5), 1, 'first');
            idx_end   = find(time_ms_eeg <= (wc + step_size_s/2 + 1e-5), 1, 'last');
            
            % Average the instantaneous power within the 100ms window per trial
            powA_all(:, w, :, b) = mean(powA_inst(:, idx_start:idx_end, :), 2, 'omitnan');
            powB_all(:, w, :, b) = mean(powB_inst(:, idx_start:idx_end, :), 2, 'omitnan');
        end
    end
end