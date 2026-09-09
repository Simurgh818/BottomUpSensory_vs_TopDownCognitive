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
window_size_ms = 100; 
step_size_ms = 100;
bg_windows = [0, 0.5]; 

bands = struct('alpha', [8 13], 'beta', [13 30]);
band_names = fieldnames(bands);
win_sizes_ms = struct('alpha', 250, 'beta', 125);

pairs = {
    {'BLT', 'P1',        [0.9, 1.5], 'Tactile vs Cued (-100 to 500ms)'},
    {'BLA', 'P1',        [0.4, 1.0], 'Auditory vs Cued (-100 to 500ms)'},
    {'P1',  'P2',        [0.9, 3.0], 'Cued vs Unpred (-100 to 2000ms)'},       
    {'P1',  'P2_500',    [0.9, 1.5], 'Cued vs Unpred 500 (-100 to 500ms)'},
    {'P1',  'P2_2000',   [2.4, 3.0], 'Cued vs Unpred 2000 (-100 to 500ms)'}, 
    {'P1',  'P3',        [0.9, 1.5], 'Cued vs Rand Cued (-100 to 500ms)'},
    {'P1',  'P3_500',    [0.9, 1.5], 'Cued vs Rand 500 (-100 to 500ms)'},
    {'P1',  'P3_missing', [0.9, 1.5], 'Cued vs Rand Missing (-100 to 500ms)'}
};

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
if isempty(set_files)
    error('No .set files found in %s', first_cond_dir);
end
names_cell = cellstr({set_files.name})';
subj_numbers = zeros(length(names_cell), 1);
for i = 1:length(names_cell)
    num_match = regexp(names_cell{i}, '\d+', 'match');
    if ~isempty(num_match)
        subj_numbers(i) = str2double(num_match{end});
    else
        subj_numbers(i) = i;
    end
end
[~, sort_idx] = sort(subj_numbers);
names_sorted = names_cell(sort_idx);
num_subjects = length(names_sorted);

get_clean_name = @(c) strrep(strrep(strrep(strrep(strrep(strrep(c, ...
    'BLA', 'Auditory'), 'BLT', 'Tactile'), 'P1', 'Cued'), ...
    'P2', 'Unpred.'), 'P3', 'Rand. Cued'), '_', ' ');

% =========================================================================
% 3.5 DATA STORAGE INITIALIZATION
% =========================================================================
GROUP_DIFF_CONN = cell(num_subjects, length(pairs)); 
GROUP_TIME_AXIS = cell(num_subjects, length(pairs)); 
GROUP_TRIAL_POWER = cell(num_subjects, length(pairs), length(band_names)); % NEW: Power Tensor

sample_EEG = pop_loadset('filename', names_sorted{1}, 'filepath', first_cond_dir, 'loadmode', 'info');
chanlocs = sample_EEG.chanlocs;
all_channels_str = {chanlocs.labels};

% --- 4. Outer Loop: Subjects (PARALLELIZED) ---
parfor target_subj = 1:num_subjects
    
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
    fprintf('Worker Processing SUBJECT %d / %d (%s)\n', target_subj, num_subjects, subj_id);
    fprintf('======================================================\n');
    
    [subject_data, time_ms_eeg, fs, ~] = load_subject_eeg(input_path, conditions, num_ch, subj_id);
    subj_dir = fullfile(output_path, subj_id);
    if ~exist(subj_dir, 'dir'), mkdir(subj_dir); end
    
    temp_subj_diff = cell(length(pairs), 1);
    temp_subj_time = cell(length(pairs), 1);
    temp_power_cell = cell(length(pairs), length(band_names));
    
    % --- 5. Middle Loop: Condition Pairs ---
    for p = 1:length(pairs)
        condA = pairs{p}{1};
        condB = pairs{p}{2};
        t_win = pairs{p}{3};      
        state_name = pairs{p}{4}; 
        
        cleanA = get_clean_name(condA);
        cleanB = get_clean_name(condB);
        
        if ~isfield(subject_data, condA) || ~isfield(subject_data, condB)
            continue;
        end
        
        trialsA_raw = subject_data.(condA);
        trialsB_raw = subject_data.(condB);
        if isempty(trialsA_raw) || isempty(trialsB_raw), continue; end
        
        % --- NEW: Extract Single-Trial Windowed Power for Group Permutation ---
        fprintf('   [%s] -> Extracting single-trial instantaneous power...\n', subj_id);
        [powA_all, powB_all, p_centers] = extract_trial_power_windows(trialsA_raw, trialsB_raw, t_win, time_ms_eeg, fs, step_size_ms, bands);
        
        for b = 1:length(band_names)
            temp_power_cell{p, b} = {powA_all(:,:,:,b), powB_all(:,:,:,b)};
        end
        
        fprintf('   [%s] -> Processing Connectivity Pair: %s vs %s...\n', subj_id, condA, condB);
        diff_state_cell = cell(1, length(band_names));
        sA_cell = cell(1, length(band_names));
        sB_cell = cell(1, length(band_names));
        win_centers_arr = [];
        
       % --- 6. Inner Loop: EEG Frequency Bands (Connectivity) ---
        for b = 1:length(band_names)
            current_band = band_names{b};
            f_range = bands.(current_band);
            
            trialsA_filt = filter_trials_band(trialsA_raw, f_range, fs);
            trialsB_filt = filter_trials_band(trialsB_raw, f_range, fs);
            
            [sA, sB, spcA, spcB, wc, k, pcl] = compute_dynamic_connectivity(trialsA_filt, trialsB_filt, ...
                t_win, time_ms_eeg, fs, window_size_ms, step_size_ms, bg_windows);
            
            diff_state_cell{b} = sA - sB; 
            sA_cell{b} = sA;
            sB_cell{b} = sB;
            if b == 1, win_centers_arr = wc; end
        end
        
        temp_subj_diff{p} = diff_state_cell;
        temp_subj_time{p} = win_centers_arr;
        
        % Within-subject topoplots can still be generated here
        plot_within_subj_topos(sA_cell, sB_cell, win_centers_arr, window_size_ms, band_names, cleanA, cleanB, state_name, chanlocs, subj_dir, subj_id);
    end
    
    GROUP_DIFF_CONN(target_subj, :) = temp_subj_diff;
    GROUP_TIME_AXIS(target_subj, :) = temp_subj_time;
    GROUP_TRIAL_POWER(target_subj, :, :) = temp_power_cell;
end
disp('All parallel subject processing complete!');

% =========================================================================
% --- 8. BETWEEN-SUBJECT GROUP LEVEL STATISTICS ---
% =========================================================================
disp('Calculating Between-Subject Group Statistics...');
group_out_dir = fullfile(output_path, 'Group_Level_Results');
if ~exist(group_out_dir, 'dir'), mkdir(group_out_dir); end

n_perms = 1000;
alpha_level = 0.05;

for p = 1:length(pairs)
    condA = pairs{p}{1}; 
    condB = pairs{p}{2};
    cleanA = get_clean_name(condA); 
    cleanB = get_clean_name(condB);
    state_name = pairs{p}{4};
    
    % Safe extraction of t_axis across 1D or 2D cell formats
    if ismatrix(GROUP_TIME_AXIS) && size(GROUP_TIME_AXIS, 1) >= 1
        t_axis = GROUP_TIME_AXIS{1, p};
    else
        t_axis = GROUP_TIME_AXIS{p};
    end
    
    if isempty(t_axis), continue; end
    
    % ---------------------------------------------------------------------
    % 1. Connectivity Networks Group Statistics (T-Test vs 0)
    % ---------------------------------------------------------------------
    for b = 1:length(band_names)
        band = band_names{b};
        valid_subjs = 0;
        group_tensor = [];
        
        for s = 1:num_subjects
            if ~isempty(GROUP_DIFF_CONN{s, p})
                valid_subjs = valid_subjs + 1;
                group_tensor(valid_subjs, :, :, :) = GROUP_DIFF_CONN{s, p}{b};
            end
        end
        
        if valid_subjs > 1
            grand_avg_net = squeeze(mean(group_tensor, 1, 'omitnan'));
            [~, p_values, ~, ~] = ttest(group_tensor, 0, 'Alpha', 0.10, 'Dim', 1);
            
            % if exist('plot_group_level_networks', 'file') == 2
            %     plot_group_level_networks(grand_avg_net, squeeze(p_values), t_axis, ...
            %         cleanA, cleanB, state_name, band, all_channels_str, chanlocs, group_out_dir);
            % end
        end
    end
    
    % ---------------------------------------------------------------------
    % 2. Group-Level Band Power Index & Beta/Alpha Ratio Index Figures
    % ---------------------------------------------------------------------
    fprintf('\nGenerating Group-Level Power Index Figures for %s vs %s...\n', cleanA, cleanB);
    alpha_idx = find(strcmpi(band_names, 'alpha'));
    beta_idx  = find(strcmpi(band_names, 'beta'));

    if ~isempty(alpha_idx) && ~isempty(beta_idx)
        pow_A_group = cell(num_subjects, 2);
        pow_B_group = cell(num_subjects, 2);

        for s = 1:num_subjects
            if ~isempty(GROUP_TRIAL_POWER{s, p, alpha_idx}) && ~isempty(GROUP_TRIAL_POWER{s, p, beta_idx})
                pow_A_group{s, 1} = GROUP_TRIAL_POWER{s, p, alpha_idx}{1}; % Alpha Cond A
                pow_B_group{s, 1} = GROUP_TRIAL_POWER{s, p, alpha_idx}{2}; % Alpha Cond B
                pow_A_group{s, 2} = GROUP_TRIAL_POWER{s, p, beta_idx}{1};  % Beta Cond A
                pow_B_group{s, 2} = GROUP_TRIAL_POWER{s, p, beta_idx}{2};  % Beta Cond B
            end
        end

        % Generates both Band Power Index (Alpha & Beta) and Beta/Alpha Ratio Index
        plot_group_power_index(pow_A_group, pow_B_group, t_axis, ...
            cleanA, cleanB, state_name, chanlocs, n_perms, alpha_level, group_out_dir);
    end
    
end
disp('Group-Level analytical extraction complete!');