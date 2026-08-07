function trials_filt = filter_trials_band(trials_raw, f_range, fs)
    % trials_raw is [Channels x Time x Trials]
    num_trials = size(trials_raw, 3);
    trials_filt = zeros(size(trials_raw));
    
    for tr = 1:num_trials
        % Transpose to [Time x Channels] for bandpass
        trial_data = trials_raw(:, :, tr)';
        try
            filt_data = bandpass(trial_data, f_range, fs)';
        catch
            warning('bandpass failed. Returning unfiltered data.');
            filt_data = trial_data';
        end
        trials_filt(:, :, tr) = filt_data;
    end
end