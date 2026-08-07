%% main_eeg_connectivity.m
clear; clc; close all;

% --- 1. Define Paths ---
if exist('I:\', 'dir')
    input_path = 'I:\My Drive\Data\New Data\EEG epoched\';
    base_output_path = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper';
elseif exist('G:\', 'dir')
    input_path = 'G:\My Drive\Data\New Data\EEG epoched\';
    base_output_path = 'C:\Users\sdabiri\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper';
else
    error('Unknown system: Cannot determine input and output paths.');
end
output_path = fullfile(base_output_path, 'functional_connectivity');

% --- 2. Configuration Parameters ---
conditions = {'BLA','BLT','P1','P2','P3'};
num_ch = 32;
n_perms = 1000;
window_size_ms = 100;
step_size_ms = 50;

% Define Condition Pairs {CondA, CondB, [Stim_Win], [Cog_Win]}
pairs = {
    {'BLA', 'P1', [0.5, 0.85], [0.85, 1.0]},
    {'BLT', 'P1', [1.0, 1.35], [1.35, 1.5]},
    {'P1', 'P2_500', [1.0, 1.35], [1.35, 1.5]},    
    {'P1', 'P2_2000', [2.5, 2.85], [2.85, 3.0]}, 
    {'P1', 'P3_500', [1.0, 1.35], [1.35, 1.5]},
    {'P1', 'P3_missing', [0.5, 0.85], [1, 1.5]}
};
state_names = {'Stimulus State', 'Cognitive State'};

% Define Canonical EEG Bands
bands = struct('delta',[1 4], 'theta',[4 8], 'alpha',[8 13], 'beta',[13 30], 'gamma',[30 50]);
band_names = fieldnames(bands);

% Start Parallel Pool
target_workers = 5; 
current_pool = gcp('nocreate');
if isempty(current_pool)
    parpool(target_workers);
elseif current_pool.NumWorkers ~= target_workers
    delete(current_pool);
    parpool(target_workers);
end

% --- 3. Determine Number of Subjects & Get Filenames ---
first_cond_dir = fullfile(input_path, conditions{1});
set_files = dir(fullfile(first_cond_dir, '*.set'));
names_sorted = sort(cellstr({set_files.name})'); %
num_subjects = length(set_files);

if num_subjects == 0
    error('No .set files found in %s', first_cond_dir);
end

% --- 4. Outer Loop: Subjects (Serial to save RAM) ---
for target_subj = 1:2 % num_subjects
    
    % DYNAMIC SUBJECT ID EXTRACTION[cite: 1]
    base_file = names_sorted{target_subj};
    % Parse exactly what is between "Avg" and ".set" (e.g., 'BOS10')
    tok = regexp(base_file, 'Avg(.*?)\.set', 'tokens');
    if ~isempty(tok)
        subj_id = tok{1}{1}; 
    else
        % Fallback: take the final string snippet if 'Avg' doesn't exist
        [~, fname, ~] = fileparts(base_file);
        tmp = split(fname, ' ');
        subj_id = tmp{end};
    end
    
    fprintf('\n======================================================\n');
    fprintf('PROCESSING SUBJECT %d / %d (%s)\n', target_subj, num_subjects, subj_id);
    fprintf('======================================================\n');
    
    % LOAD ONLY THE CURRENT SUBJECT INTO MEMORY
    [subject_data, time_ms_eeg, fs, all_channels_str] = load_subject_eeg(input_path, conditions, num_ch, target_subj);
    
    % Use the dynamic subj_id for the output folder
    subj_dir = fullfile(output_path, subj_id);
    if ~exist(subj_dir, 'dir'), mkdir(subj_dir); end
    
    % --- 5. Middle Loop: Condition Pairs ---
    for p = 1:length(pairs)
        condA = pairs{p}{1};
        condB = pairs{p}{2};
        
        if ~isfield(subject_data, condA) || ~isfield(subject_data, condB)
            continue;
        end
        
        % Extract raw trials for this specific subject directly from the struct
        trialsA_raw = subject_data.(condA);
        trialsB_raw = subject_data.(condB);
        
        if isempty(trialsA_raw) || isempty(trialsB_raw), continue; end

        fprintf('   -> Generating 5-Band Power Bar Charts...\n');
        
        % --- NEW: Generate Multi-Band Power Bar Charts (5-Row Figure) ---
        for state_idx = 1:2
            t_win_power = pairs{p}{2 + state_idx}; 
            active_state_power = state_names{state_idx};
            
            plot_band_power_bars(trialsA_raw, trialsB_raw, t_win_power, time_ms_eeg, fs, ...
                condA, condB, active_state_power, bands, all_channels_str, subj_dir, subj_id);
        end
        
        fprintf('   -> Deploying Pair: %s vs %s to Parallel Workers...\n', condA, condB);
        
        % --- 6. Inner PARALLEL Loop: EEG Frequency Bands ---
        parfor b = 1:length(band_names)
            current_band = band_names{b};
            f_range = bands.(current_band);
            
            % Create band directory
            band_dir = fullfile(subj_dir, current_band);
            if ~exist(band_dir, 'dir'), mkdir(band_dir); end
            
            % Filter Trials into the current frequency band
            trialsA_filt = filter_trials_band(trialsA_raw, f_range, fs);
            trialsB_filt = filter_trials_band(trialsB_raw, f_range, fs);
            
            % Loop over Stimulus and Cognitive States
            for state_idx = 1:2
                t_win = pairs{p}{2 + state_idx}; 
                active_state = state_names{state_idx};
                
                % Compute Connectivity (dPCA + Sliding Correlation + Fisher Z Subtraction)
                [sub_ch_A, sub_ch_B, sub_pc_A, sub_pc_B, window_centers, k_opt, pc_labels] = ...
                    compute_dynamic_connectivity(trialsA_filt, trialsB_filt, t_win, time_ms_eeg, fs, window_size_ms, step_size_ms, n_perms);
                
                % Plot and Save Figures (Now passing 'subj_id' forward)
                plot_dynamic_networks(sub_ch_A, sub_ch_B, sub_pc_A, sub_pc_B, window_centers, ...
                    trialsA_filt, trialsB_filt, time_ms_eeg, t_win, fs, window_size_ms, step_size_ms, ...
                    condA, condB, active_state, current_band, k_opt, pc_labels, all_channels_str, band_dir, subj_id);
            end
        end
    end
    
    % --- 7. CLEAR MEMORY ---
    clear subject_data trialsA_raw trialsB_raw;
end

disp('All parallel processing complete!');