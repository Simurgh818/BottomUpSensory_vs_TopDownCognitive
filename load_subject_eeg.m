function [subject_data, time_ms_eeg, fs, all_channels_str] = load_subject_eeg(input_path, conditions, num_ch, subj_id)
    % Initialize outputs
    subject_data = struct();
    time_ms_eeg = [];
    fs = [];
    all_channels_str = arrayfun(@(x) sprintf('Ch%d', x), 1:num_ch, 'UniformOutput', false);
    
    for c = 1:length(conditions)
        condition = conditions{c};
        in_dir = fullfile(input_path, condition);
        
        % --- NEW: Explicitly search for the file containing the subj_id ---
        search_pattern = sprintf('*%s.set', subj_id);
        set_files = dir(fullfile(in_dir, search_pattern));
        
        if isempty(set_files)
            warning('Missing data for %s in condition %s. Skipping.', subj_id, condition);
            continue;
        end
        
        % Take the exact matched file
        file_name = set_files(1).name;
        fprintf('   -> Loading %s for Condition: %s\n', file_name, condition);
        
        % Load the file
        EEG = pop_loadset(fullfile(in_dir, file_name));
        
        % Extract sampling rate on the first load
        if isempty(fs)
            fs = EEG.srate; 
        end
        
        % Apply the odd/even indexing logic pipeline
        if ismember(condition, {'BLA','P1','P2','P3'})
            epoch_trials = 1:2:EEG.trials;
        else
            epoch_trials = 2:2:EEG.trials;
        end
        
        % Extract trials [Channels x Time x Trials]
        trials = EEG.data(1:num_ch, :, epoch_trials);
        subject_data.(condition) = trials;
        
        % Build time vector on the first load
        if isempty(time_ms_eeg)
            nTime = size(trials, 2);
            time_ms_eeg = linspace(0, 3.5, nTime); % Based on 0 to 3.5 sec epochs
        end
    end
end