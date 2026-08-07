function [subject_data, time_ms_eeg, fs, all_channels_str] = load_subject_eeg(input_path, conditions, num_ch, target_subj)
    % This function loads only ONE subject into memory to prevent RAM exhaustion.
    subject_data = struct();
    time_ms_eeg = [];
    fs = [];
    all_channels_str = {};
    
    for c = 1:length(conditions)
        condition = conditions{c};
        in_dir = fullfile(input_path, condition);
        set_files = dir(fullfile(in_dir, '*.set'));
        names_sorted = sort(cellstr({set_files.name})');
        
        % Check if the subject exists for this condition
        if isempty(names_sorted) || target_subj > length(names_sorted)
            continue; 
        end
        
        % Excel specific loading logic (P2 and P3)
        if strcmp(condition, 'P2')
            epoch_trials_p2_500ms = readmatrix(fullfile(input_path, 'Indexes for P2.xlsx'), 'Sheet', 'Audio onset with 500 ms tactile');
            epoch_trials_p2_2000ms = readmatrix(fullfile(input_path, 'Indexes for P2.xlsx'), 'Sheet', 'Audio onset with 2000 ms tactil');
        elseif strcmp(condition, 'P3')
            epoch_trials_p3_500ms = readmatrix(fullfile(input_path, 'Indexes for P3.xlsx'), 'Sheet', 'Audio onset with 500 ms tactile');
            epoch_trials_p3_missing = readmatrix(fullfile(input_path, 'Indexes for P3.xlsx'), 'Sheet', 'Audio onset with missing tactil');
        end
        
        % Load ONLY the target subject
        file_to_load = names_sorted{target_subj}; 
        fprintf('   -> Loading %s for Condition: %s\n', file_to_load, condition);
        EEG = pop_loadset('filename', file_to_load, 'filepath', in_dir);
        
        if isempty(fs), fs = EEG.srate; end
        if isempty(time_ms_eeg), time_ms_eeg = linspace(0, 3.5, size(EEG.data, 2)); end
        if isempty(all_channels_str)
            if isfield(EEG, 'chanlocs'), all_channels_str = {EEG.chanlocs.labels};
            else, all_channels_str = arrayfun(@(x) sprintf('Ch%d', x), 1:num_ch, 'UniformOutput', false); end
        end
        
        % Store directly as a 3D matrix (No cell array nesting needed!)
        if ismember(condition, {'BLA', 'BLT', 'P1'})
            subject_data.(condition) = EEG.data(:, :, 1:2:EEG.trials);
        elseif strcmp(condition, 'P2')
            subject_data.(condition) = EEG.data(:, :, 1:2:EEG.trials);
            subject_data.P2_500      = EEG.data(:, :, epoch_trials_p2_500ms);
            subject_data.P2_2000     = EEG.data(:, :, epoch_trials_p2_2000ms);
        elseif strcmp(condition, 'P3')
            subject_data.(condition)   = EEG.data(:, :, 1:2:EEG.trials);
            subject_data.P3_500        = EEG.data(:, :, epoch_trials_p3_500ms);
            subject_data.P3_missing    = EEG.data(:, :, epoch_trials_p3_missing);
        end
    end
end