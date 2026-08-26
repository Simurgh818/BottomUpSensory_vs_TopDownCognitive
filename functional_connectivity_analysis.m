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

% UPDATED: 100ms window size and 100ms steps to create continuous non-overlapping frames
window_size_ms = 100; 
step_size_ms = 100;
bg_windows = [0, 0.5]; % Universal Resting State Baseline Windows (0 to 0.5s)

% Limit to just Alpha and Beta
bands = struct('alpha', [8 13], 'beta', [13 30]);
band_names = fieldnames(bands);

% Set FFT window sizes to 100ms to match the topological windowing
win_sizes_ms = struct('alpha', 250, 'beta', 125);

% --- NEW: Single Continuous Evaluation Window (-100ms to +500ms) ---
% Format: {CondA, CondB, [Eval-Win], Display-Name}
% Note: Tactile/Cued stimuli hit at 1.0s, Auditory hits at 0.5s, P2_2000 hits at 2.5s.
pairs = {
    {'BLT', 'P1',        [0.9, 1.5], 'Tactile vs Cued (-100 to 500ms)'},
    {'BLA', 'P1',        [0.4, 1.0], 'Auditory vs Cued (-100 to 500ms)'},
    {'P1',  'P2',        [0.9, 3.0], 'Cued vs Unpred (-100 to 2000ms)'},       % 21 frames (spanning 0.9s to 3.0s)
    {'P1',  'P2_500',    [0.9, 1.5], 'Cued vs Unpred 500 (-100 to 500ms)'},
    {'P1',  'P2_2000',   [2.4, 3.0], 'Cued vs Unpred 2000 (-100 to 500ms)'}, 
    {'P1',  'P3',        [0.9, 1.5], 'Cued vs Rand Cued (-100 to 500ms)'},
    {'P1',  'P3_500',    [0.9, 1.5], 'Cued vs Rand 500 (-100 to 500ms)'},
    {'P1',  'P3_missing', [0.9, 1.5], 'Cued vs Rand Missing (-100 to 500ms)'}
};

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
if isempty(set_files)
    error('No .set files found in %s', first_cond_dir);
end
names_cell = cellstr({set_files.name})';

% Extract numerical subject index
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
% num_subjects = 1; % <--- Uncomment to test just 1 subject

get_clean_name = @(c) strrep(strrep(strrep(strrep(strrep(strrep(c, ...
    'BLA', 'Auditory'), 'BLT', 'Tactile'), 'P1', 'Cued'), ...
    'P2', 'Unpred.'), 'P3', 'Rand. Cued'), '_', ' ');

% =========================================================================
% 3.5 DATA STORAGE & SPATIAL INITIALIZATION
% =========================================================================
% Structure simplified since we only have 1 continuous window sequence per pair
GROUP_DIFF_CONN = cell(num_subjects, length(pairs)); 
GROUP_TIME_AXIS = cell(num_subjects, length(pairs)); 

% Extract True 10-20 Channel Locations for the Topoplots early
sample_EEG = pop_loadset('filename', names_sorted{1}, 'filepath', first_cond_dir, 'loadmode', 'info');
chanlocs = sample_EEG.chanlocs;
all_channels_str = {chanlocs.labels};

% --- 4. Outer Loop: Subjects (PARALLELIZED) ---
parfor target_subj = 1:1 %num_subjects
    
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
    fprintf('Worker Processing SUBJECT %d / %d (%s)\n', target_subj, num_subjects, subj_id);
    fprintf('======================================================\n');
    
    [subject_data, time_ms_eeg, fs, ~] = load_subject_eeg(input_path, conditions, num_ch, subj_id);
    
    subj_dir = fullfile(output_path, subj_id);
    if ~exist(subj_dir, 'dir'), mkdir(subj_dir); end
    
    temp_subj_diff = cell(length(pairs), 1);
    temp_subj_time = cell(length(pairs), 1);
    
% =========================================================================
    % 4.5. BUILD THE GLOBAL BROADBAND MANIFOLD (ALL CONDITIONS)
    % =========================================================================
    fprintf('   [%s] -> Building Global dPCA Manifold...\n', subj_id);
    
    % Identify all available conditions for this subject
    avail_conds = fieldnames(subject_data);
    stack_avg = [];
    avg_broadband = struct();
    
    for c = 1:length(avail_conds)
        c_name = avail_conds{c};
        trials = subject_data.(c_name);
        if ~isempty(trials) && size(trials, 3) > 0
            % Trial Average
            avg = mean(trials, 3, 'omitnan');
            avg_broadband.(c_name) = avg;
            
            % Stack for global manifold [Channels x Time x Conditions]
            stack_avg = cat(3, stack_avg, avg);
        end
    end
    
    % Run MAP test and SVD/dPCA fallback
    data_2d = reshape(stack_avg, num_ch, []);

    % Pass 'true' to suppress plotting inside the parfor loop!
    [k_opt, ~] = velicer_map((data_2d - mean(data_2d, 2))', true);
    if k_opt < 2, k_opt = 2; end
    fprintf('      -> Velicer MAP assigned %d dimensions for Global Manifold.\n', k_opt);
    
    try
        X_dpca = bsxfun(@minus, stack_avg, mean(stack_avg(:,:), 2));
        [W_dpca, ~, ~] = dpca(X_dpca, k_opt);
        W = W_dpca'; % [k x Channels]
    catch
        fprintf('      -> dPCA failed or missing. Using SVD fallback.\n');
        [~, ~, V] = svd((data_2d - mean(data_2d, 2))', 'econ');
        W = V(:, 1:k_opt)';
    end

    % --- 5. Middle Loop: Condition Pairs ---
    for p = 1:length(pairs)
        condA = pairs{p}{1};
        condB = pairs{p}{2};
        t_win = pairs{p}{3};      % Eval Window
        state_name = pairs{p}{4}; % Display Name
        
        cleanA = get_clean_name(condA);
        cleanB = get_clean_name(condB);
        
        if ~isfield(subject_data, condA) || ~isfield(subject_data, condB)
            continue;
        end
        
        trialsA_raw = subject_data.(condA);
        trialsB_raw = subject_data.(condB);
        if isempty(trialsA_raw) || isempty(trialsB_raw), continue; end
        
        fprintf('   [%s] -> Generating Continuous Band Power Topoplots...\n', subj_id);
        
        plot_band_power_topos(trialsA_raw, trialsB_raw, t_win, time_ms_eeg, fs, ...
            win_sizes_ms, step_size_ms, condA, condB, state_name, bands, chanlocs, subj_dir, subj_id);
        
        fprintf('   [%s] -> Processing Pair: %s vs %s across Alpha & Beta...\n', subj_id, condA, condB);
        
        % --- NEW: GENERATE dPCA TRAJECTORY & SUBSPACE FILMSTRIPS ---
        fprintf('   [%s] -> Generating dPCA Trajectory & Subspace Filmstrips...\n', subj_id);
        
        avgA = avg_broadband.(condA);
        avgB = avg_broadband.(condB);
        
        plot_dpca_subspace_filmstrip(avgA, avgB, W, k_opt, t_win, time_ms_eeg, ...
            window_size_ms, step_size_ms, bg_windows, cleanA, cleanB, ...
            all_channels_str, subj_dir, subj_id);

        diff_state_cell  = cell(1, length(band_names));
        sA_cell          = cell(1, length(band_names));
        sB_cell          = cell(1, length(band_names));
        win_centers_arr  = [];
        
       % --- 6. Inner Loop: EEG Frequency Bands ---
        for b = 1:length(band_names)
            current_band = band_names{b};
            f_range = bands.(current_band);
            
            band_dir = fullfile(subj_dir, current_band);
            if ~exist(band_dir, 'dir'), mkdir(band_dir); end
            
            trialsA_filt = filter_trials_band(trialsA_raw, f_range, fs);
            trialsB_filt = filter_trials_band(trialsB_raw, f_range, fs);
            
            % --- GENERATE 3-ROW CORRELATION FILMSTRIPS (Notebook Implementation) ---
            plot_correlation_filmstrip(trialsA_filt, t_win, time_ms_eeg, window_size_ms, step_size_ms, ...
                bg_windows, cleanA, current_band, all_channels_str, band_dir, subj_id);
            plot_correlation_filmstrip(trialsB_filt, t_win, time_ms_eeg, window_size_ms, step_size_ms, ...
                bg_windows, cleanB, current_band, all_channels_str, band_dir, subj_id);
            
            % --- Dimensionality Reduction & Network Dynamics ---
            [sA, sB, spcA, spcB, wc, k, pcl] = compute_dynamic_connectivity(trialsA_filt, trialsB_filt, ...
                t_win, time_ms_eeg, fs, window_size_ms, step_size_ms, bg_windows);
            
            diff_state_cell{b} = sA - sB; 
            sA_cell{b} = sA;
            sB_cell{b} = sB;
            if b == 1, win_centers_arr = wc; end
        end
        
        temp_subj_diff{p} = diff_state_cell;
        temp_subj_time{p} = win_centers_arr;
       
        % Generate the Within-Subject Connectivity Topoplots (Delta r Index)
        % --- UPDATED: Added window_size_ms to the function call ---
        plot_within_subj_topos(sA_cell, sB_cell, win_centers_arr, window_size_ms, band_names, cleanA, cleanB, state_name, chanlocs, subj_dir, subj_id);
    end
    
    GROUP_DIFF_CONN(target_subj, :) = temp_subj_diff;
    GROUP_TIME_AXIS(target_subj, :) = temp_subj_time;
end
disp('All parallel subject processing complete!');

% --- 8. BETWEEN-SUBJECT GROUP LEVEL STATISTICS ---
disp('Calculating Between-Subject Group Statistics...');
group_out_dir = fullfile(output_path, 'Group_Level_Results');
if ~exist(group_out_dir, 'dir'), mkdir(group_out_dir); end

for p = 1:length(pairs)
    condA = pairs{p}{1}; condB = pairs{p}{2};
    cleanA = get_clean_name(condA); cleanB = get_clean_name(condB);
    state_name = pairs{p}{4};
    
    t_axis = GROUP_TIME_AXIS{1, p};
    if isempty(t_axis), continue; end
    
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
            
            plot_group_level_networks(grand_avg_net, squeeze(p_values), t_axis, cleanA, cleanB, state_name, band, all_channels_str, chanlocs, group_out_dir);
        end
    end
end
disp('Group-Level extraction complete!');