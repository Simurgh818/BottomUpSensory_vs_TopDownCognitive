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

% --- 3. Determine Number of Subjects & Get Filenames (Natural Numerical Order) ---
first_cond_dir = fullfile(input_path, conditions{1});
set_files = dir(fullfile(first_cond_dir, '*.set'));
if isempty(set_files)
    error('No .set files found in %s', first_cond_dir);
end

names_cell = cellstr({set_files.name})';

% Extract numerical subject index (e.g. BOS2 -> 2, BOS10 -> 10)
subj_numbers = zeros(length(names_cell), 1);
for i = 1:length(names_cell)
    num_match = regexp(names_cell{i}, '\d+', 'match');
    if ~isempty(num_match)
        subj_numbers(i) = str2double(num_match{end}); % Grabs trailing digits
    else
        subj_numbers(i) = i;
    end
end

[~, sort_idx] = sort(subj_numbers);
names_sorted = names_cell(sort_idx);
num_subjects = length(names_sorted);

if num_subjects == 0
    error('No .set files found in %s', first_cond_dir);
end

get_clean_name = @(c) strrep(strrep(strrep(strrep(strrep(strrep(c, ...
    'BLA', 'Auditory'), 'BLT', 'Tactile'), 'P1', 'Cued'), ...
    'P2', 'Unpred.'), 'P3', 'Rand. Cued'), '_', ' ');

% =========================================================================
% NEW: 3.5 DATA STORAGE ARCHITECTURE INITIALIZATION
% This explicitly sets up the group-level tensors before the subject loop
% Dimensions: {Subject, Pair, State} -> cell array of bands [Ch x Ch x Win]
% =========================================================================
GROUP_DIFF_CONN = cell(num_subjects, length(pairs), 2); 
GROUP_TIME_AXIS = cell(length(pairs), 2);


% --- 4. Outer Loop: Subjects (Serial to save RAM) ---
for target_subj = 1:num_subjects
    
    % DYNAMIC SUBJECT ID EXTRACTION
    base_file = names_sorted{target_subj};
    tok = regexp(base_file, 'Avg(.*?)\.set', 'tokens');
    if ~isempty(tok)
        subj_id = tok{1}{1}; 
    else
        [~, fname, ~] = fileparts(base_file);
        tmp = split(fname, ' ');
        subj_id = tmp{end};
    end
    
    fprintf('\n======================================================\n');
    fprintf('PROCESSING SUBJECT %d / %d (%s)\n', target_subj, num_subjects, subj_id);
    fprintf('======================================================\n');

    % LOAD ONLY THE CURRENT SUBJECT INTO MEMORY
    [subject_data, time_ms_eeg, fs, all_channels_str] = load_subject_eeg(input_path, conditions, num_ch, subj_id);
    
    subj_dir = fullfile(output_path, subj_id);
    if ~exist(subj_dir, 'dir'), mkdir(subj_dir); end
    
    % --- 5. Middle Loop: Condition Pairs ---
    for p = 1:length(pairs)
        condA = pairs{p}{1};
        condB = pairs{p}{2};
        cleanA = get_clean_name(condA);
        cleanB = get_clean_name(condB);
        
        if ~isfield(subject_data, condA) || ~isfield(subject_data, condB)
            continue;
        end
        
        trialsA_raw = subject_data.(condA);
        trialsB_raw = subject_data.(condB);
        
        if isempty(trialsA_raw) || isempty(trialsB_raw), continue; end
        
        fprintf('   -> Generating 5-Band Power Bar Charts...\n');
        for state_idx = 1:2
            t_win_power = pairs{p}{2 + state_idx}; 
            plot_band_power_bars(trialsA_raw, trialsB_raw, t_win_power, time_ms_eeg, fs, ...
                condA, condB, state_names{state_idx}, bands, all_channels_str, subj_dir, subj_id);
        end
        
        fprintf('   -> Deploying Pair: %s vs %s to Parallel Workers...\n', condA, condB);
        
        % NEW: Preallocate containers to safely catch parfor outputs
        diff_state1_cell = cell(1, length(band_names));
        diff_state2_cell = cell(1, length(band_names));
        win_centers1_cell = cell(1, length(band_names));
        win_centers2_cell = cell(1, length(band_names));
        
        % --- 6. Inner PARALLEL Loop: EEG Frequency Bands ---
        parfor b = 1:length(band_names)
            current_band = band_names{b};
            f_range = bands.(current_band);
            
            band_dir = fullfile(subj_dir, current_band);
            if ~exist(band_dir, 'dir'), mkdir(band_dir); end
            
            trialsA_filt = filter_trials_band(trialsA_raw, f_range, fs);
            trialsB_filt = filter_trials_band(trialsB_raw, f_range, fs);
            
            % --- State 1 (Stimulus) ---
            t_win1 = pairs{p}{3}; 
            [sA1, sB1, spcA1, spcB1, wc1, k1, pcl1] = compute_dynamic_connectivity(trialsA_filt, trialsB_filt, t_win1, time_ms_eeg, fs, window_size_ms, step_size_ms, n_perms);
            plot_dynamic_networks(sA1, sB1, spcA1, spcB1, wc1, trialsA_filt, trialsB_filt, time_ms_eeg, t_win1, fs, window_size_ms, step_size_ms, condA, condB, state_names{1}, current_band, k1, pcl1, all_channels_str, band_dir, subj_id);
            
            % Save subtracted network (Test - Control) for State 1
            diff_state1_cell{b} = sB1 - sA1; 
            win_centers1_cell{b} = wc1;
            
            % --- State 2 (Cognitive) ---
            t_win2 = pairs{p}{4}; 
            [sA2, sB2, spcA2, spcB2, wc2, k2, pcl2] = compute_dynamic_connectivity(trialsA_filt, trialsB_filt, t_win2, time_ms_eeg, fs, window_size_ms, step_size_ms, n_perms);
            plot_dynamic_networks(sA2, sB2, spcA2, spcB2, wc2, trialsA_filt, trialsB_filt, time_ms_eeg, t_win2, fs, window_size_ms, step_size_ms, condA, condB, state_names{2}, current_band, k2, pcl2, all_channels_str, band_dir, subj_id);
            
            % Save subtracted network (Test - Control) for State 2
            diff_state2_cell{b} = sB2 - sA2;
            win_centers2_cell{b} = wc2;
        end
        
       % =========================================================================
        % NEW: 6.5 POST-PARFOR AGGREGATION & WITHIN-SUBJECT PROFILE
        % =========================================================================
        % Store the raw network matrices for Group-Level tests
        GROUP_DIFF_CONN{target_subj, p, 1} = diff_state1_cell;
        GROUP_DIFF_CONN{target_subj, p, 2} = diff_state2_cell;
        GROUP_TIME_AXIS{p, 1} = win_centers1_cell{1};
        GROUP_TIME_AXIS{p, 2} = win_centers2_cell{1};
        
        % Calculate the global mean connectivity difference across all bands
        mean_diff1 = zeros(length(band_names), length(win_centers1_cell{1}));
        mean_diff2 = zeros(length(band_names), length(win_centers2_cell{1}));
        
        for b = 1:length(band_names)
            % --- Aggregate State 1 (Stimulus) ---
            mat1 = abs(diff_state1_cell{b}); 
            for w = 1:size(mat1, 3)
                % Zero the diagonal to prevent self-correlation from inflating the mean
                t1 = mat1(:,:,w); t1(1:num_ch+1:end) = nan; 
                mean_diff1(b, w) = mean(t1, 'all', 'omitnan');
            end
            
            % --- Aggregate State 2 (Cognitive) ---
            mat2 = abs(diff_state2_cell{b});
            for w = 1:size(mat2, 3)
                % Zero the diagonal to prevent self-correlation from inflating the mean
                t2 = mat2(:,:,w); t2(1:num_ch+1:end) = nan;
                mean_diff2(b, w) = mean(t2, 'all', 'omitnan');
            end
        end
        
        % Generate the Within-Subject Band Profile Figures
        plot_within_subj_bands(mean_diff1, win_centers1_cell{1}, band_names, cleanA, cleanB, state_names{1}, subj_dir, subj_id);
        plot_within_subj_bands(mean_diff2, win_centers2_cell{1}, band_names, cleanA, cleanB, state_names{2}, subj_dir, subj_id);
    end
    
    % --- 7. CLEAR MEMORY ---
    clear subject_data trialsA_raw trialsB_raw;
end
disp('All parallel subject processing complete!');

% =========================================================================
% NEW: 8. BETWEEN-SUBJECT GROUP LEVEL STATISTICS
% This runs after ALL subjects have finished processing.
% =========================================================================
disp('Calculating Between-Subject Group Statistics...');
group_out_dir = fullfile(output_path, 'Group_Level_Results');
if ~exist(group_out_dir, 'dir'), mkdir(group_out_dir); end

% Extract Grand Average and T-Test per Band & State
for p = 1:length(pairs)
    condA = pairs{p}{1}; condB = pairs{p}{2};
    cleanA = get_clean_name(condA); cleanB = get_clean_name(condB);
    
    for state_idx = 1:2
        state_name = state_names{state_idx};
        t_axis = GROUP_TIME_AXIS{p, state_idx};
        if isempty(t_axis), continue; end
        
        for b = 1:length(band_names)
            band = band_names{b};
            
            % Collect all subjects into a 4D array: [Subjects x Channels x Channels x Windows]
            % Exclude empty subjects (e.g., if one subj missed a condition)
            valid_subjs = 0;
            group_tensor = [];
            for s = 1:num_subjects
                if ~isempty(GROUP_DIFF_CONN{s, p, state_idx})
                    valid_subjs = valid_subjs + 1;
                    group_tensor(valid_subjs, :, :, :) = GROUP_DIFF_CONN{s, p, state_idx}{b};
                end
            end
            
            if valid_subjs > 1
                % Calculate Grand Average and Mass-Univariate T-Test against 0
                grand_avg_net = squeeze(mean(group_tensor, 1, 'omitnan'));
                [~, p_values, ~, stats] = ttest(group_tensor, 0, 'Alpha', 0.05, 'Dim', 1);
                t_scores = squeeze(stats.tstat);
                p_values = squeeze(p_values);
                
                % Optional: Pass `grand_avg_net` and `p_values` into a group plotting function
                % plot_group_level_networks(grand_avg_net, p_values, t_axis, cleanA, cleanB, state_name, band, all_channels_str, group_out_dir);
            end
        end
    end
end
disp('Group-Level extraction complete!');