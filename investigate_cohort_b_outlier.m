%% =========================================================================
%% OUTLIER HUNTER & TRIAL-LEVEL SCREENING (SUBJECT 3 & SUBJECT 14)
%% =========================================================================
clear; clc; close all;

% --- 1. SETUP PATHS ---
if exist('I:\', 'dir')
    input_path = 'I:\My Drive\Data\New Data\EEG epoched\';
    base_output_path = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper\';
elseif exist('H:\', 'dir')
    input_path = 'H:\My Drive\Data\New Data\EEG epoched\';
    base_output_path = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper\';
elseif exist('G:\', 'dir')
    input_path = 'G:\My Drive\Data\New Data\EEG epoched\';
    base_output_path = 'C:\Users\sdabiri\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper\';
else
    error('Unknown system paths.');
end

raw_results_file = fullfile(base_output_path, 'ch-ch_cov', 'group', 'Raw', 'P1', 'results.mat');
if ~exist(raw_results_file, 'file')
    error('Could not find results.mat at: %s.', raw_results_file);
end

% --- 2. LOAD PRECOMPUTED GAINS ---
load(raw_results_file, 'G_splits'); % [nCond x nRep x num_subjects]
conditions = {'P1', 'P2_500', 'P3_500', 'P2_2000', 'P3_missing', 'BLA'};
idx_p3_500 = find(strcmp(conditions, 'P3_500'));
idx_p3_mis = find(strcmp(conditions, 'P3_missing'));
[nCond, nRep, num_subjs] = size(G_splits);
subj_G_matrix = squeeze(mean(G_splits, 2, 'omitnan'))'; % [num_subjs x nCond]

cohort_A_idx = intersect(4:11, 1:num_subjs);
cohort_B_idx = setdiff(1:num_subjs, cohort_A_idx);

target_subjs = [3, 14];
save_trial_cleaned_sets = true; % Set true to export cleaned .set files

%% =========================================================================
%% 3. IDENTIFY OUTLIERS ACROSS ALL SUBJECTS
%% =========================================================================
valid_subjs = [cohort_A_idx, cohort_B_idx];
excluded_subjs = setdiff(1:num_subjs, valid_subjs);
X_features = [subj_G_matrix(:, idx_p3_500), subj_G_matrix(:, idx_p3_mis)];

dist_from_median = sqrt((X_features(:,1) - median(X_features(valid_subjs, 1))).^2 + ...
                        (X_features(:,2) - median(X_features(valid_subjs, 2))).^2);
dist_from_median(excluded_subjs) = -1;
[~, sort_idx] = sort(dist_from_median, 'descend');
outlier_subj = sort_idx(1);

fprintf('========================================================================\n');
fprintf('OUTLIER DISCOVERY REPORT (ALL 14 SUBJECTS)\n');
fprintf('========================================================================\n');
for i = 1:3
    s_id = sort_idx(i);
    grp = 'B'; if ismember(s_id, cohort_A_idx), grp = 'A'; end
    fprintf('  #%d: Subject %02d (Cohort %s) | Rand 500 Gain: %6.2f | Rand Null Gain: %6.2f\n', ...
        i, s_id, grp, subj_G_matrix(s_id, idx_p3_500), subj_G_matrix(s_id, idx_p3_mis));
end
fprintf('\n=> Executing Trial-Level Screening on Target Subjects: %s\n', mat2str(target_subjs));

%% =========================================================================
%% 4. PLOT PARTICIPANT GAIN PROFILES
%% =========================================================================
fig_diag = figure('Position', [100, 100, 1600, 600], 'Name', 'Participant Gain Profiles');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
bar_cols = repmat([0.2 0.4 0.8], num_subjs, 1);
bar_cols(cohort_A_idx, :) = repmat([0.5 0.5 0.5], length(cohort_A_idx), 1);
bar_cols(target_subjs, :) = repmat([0.85 0.2 0.2], length(target_subjs), 1);

nexttile; hold on; set(gca, 'FontSize', 16);
yline(1.0, 'k--', 'LineWidth', 2);
for s = 1:num_subjs, b = bar(s, subj_G_matrix(s, idx_p3_500), 0.7); b.FaceColor = bar_cols(s, :); end
xticks(1:num_subjs); xlabel('Participant ID'); ylabel('Global Gain (G)');
title('Rand 500 (P3) Gain per Participant', 'FontSize', 18, 'FontWeight', 'bold'); grid on;

nexttile; hold on; set(gca, 'FontSize', 16);
yline(1.0, 'k--', 'LineWidth', 2);
for s = 1:num_subjs, b = bar(s, subj_G_matrix(s, idx_p3_mis), 0.7); b.FaceColor = bar_cols(s, :); end
xticks(1:num_subjs); xlabel('Participant ID'); ylabel('Global Gain (G)');
title('Rand Null (P3 Missing) Gain per Participant', 'FontSize', 18, 'FontWeight', 'bold'); grid on;
sgtitle('Gain Decomposition: Gray = Cohort A, Blue = Cohort B, Red = Screened Targets (S3 & S14)', ...
    'FontSize', 20, 'FontWeight', 'bold');

%% =========================================================================
%% 5. LOAD RAW DATA (EXACT EXCEL INDEX MATCHING)
%% =========================================================================
fprintf('\nLoading EEGLAB data using Excel Trial Indices...\n');
data_all_conds = struct();
data_all_conds.Raw = struct();
EEG_templates = struct(); 
conds_to_load = {'BLA', 'BLT', 'P1', 'P2', 'P3'};

% Pre-read trial indices from Excel
excel_p2 = fullfile(input_path, 'Indexes for P2.xlsx');
if exist(excel_p2, 'file')
    epoch_trials_p2_500  = readmatrix(excel_p2, 'Sheet', 'Audio onset with 500 ms tactile');
    epoch_trials_p2_2000 = readmatrix(excel_p2, 'Sheet', 'Audio onset with 2000 ms tactil');
else
    error('Could not find: %s', excel_p2);
end
excel_p3 = fullfile(input_path, 'Indexes for P3.xlsx');
if exist(excel_p3, 'file')
    epoch_trials_p3_500 = readmatrix(excel_p3, 'Sheet', 'Audio onset with 500 ms tactile');
    epoch_trials_p3_mis = readmatrix(excel_p3, 'Sheet', 'Audio onset with missing tactil');
else
    error('Could not find: %s', excel_p3);
end

% Remove any NaN values if read from trailing rows
epoch_trials_p2_500  = epoch_trials_p2_500(~isnan(epoch_trials_p2_500));
epoch_trials_p2_2000 = epoch_trials_p2_2000(~isnan(epoch_trials_p2_2000));
epoch_trials_p3_500  = epoch_trials_p3_500(~isnan(epoch_trials_p3_500));
epoch_trials_p3_mis  = epoch_trials_p3_mis(~isnan(epoch_trials_p3_mis));

for c_idx = 1:length(conds_to_load)
    cond_name = conds_to_load{c_idx};
    in_dir = fullfile(input_path, cond_name);
    set_files = dir(fullfile(in_dir, '*.set'));
    names_sorted = sort(cellstr({set_files.name})');
    
    if isempty(names_sorted)
        fprintf('No .set files for %s — skipping\n', cond_name);
        continue;
    end
    
    for s = 1:num_subjs
        file_to_load = names_sorted{s};
        
        % --- FORCED RAW LOADING ---
        % We deliberately ignore the _cleaned folders here so we can perform
        % trial-level screening on the pristine raw data.
        fprintf('Loading RAW Subj %d/%d (%s) for condition: %s\n', ...
            s, num_subjs, file_to_load, cond_name);
        EEG = pop_loadset('filename', file_to_load, 'filepath', in_dir);
        
        if ismember(s, target_subjs)
            EEG_templates.(cond_name){s} = EEG;
            EEG_templates.(cond_name){s}.filename_orig = file_to_load;
        end
        
        % Correct trial index assignments matching your experimental structure
        if ismember(cond_name, {'BLA', 'BLT', 'P1'})
            epoch_trials = 1:2:EEG.trials; % 60 trials for P1
            data_all_conds.Raw.(cond_name){s} = double(EEG.data(:, :, epoch_trials));
            
        elseif strcmp(cond_name, 'P2')
            data_all_conds.Raw.P2{s}      = double(EEG.data(:, :, 1:2:EEG.trials));
            data_all_conds.Raw.P2_500{s}  = double(EEG.data(:, :, epoch_trials_p2_500));  % 30 trials
            data_all_conds.Raw.P2_2000{s} = double(EEG.data(:, :, epoch_trials_p2_2000)); % 30 trials
            
        elseif strcmp(cond_name, 'P3')
            data_all_conds.Raw.P3{s}         = double(EEG.data(:, :, 1:2:EEG.trials));
            data_all_conds.Raw.P3_500{s}     = double(EEG.data(:, :, epoch_trials_p3_500)); % 30 trials (Rand 500)
            data_all_conds.Raw.P3_missing{s} = double(EEG.data(:, :, epoch_trials_p3_mis)); % 30 trials (Rand Null)
        end
    end
end

% Set electrode C3 index
chan_idx = 14; 
if exist('EEG', 'var') && isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs)
    c3_match = find(strcmpi({EEG.chanlocs.labels}, 'C3'), 1);
    if ~isempty(c3_match), chan_idx = c3_match; end
end
fprintf('\nUsing Channel Index %d (C3) for ERP waveform verification.\n', chan_idx);

nPnts  = size(data_all_conds.Raw.P1{1}, 2);
time_s = -1.000 + (0:(nPnts - 1)) / 500;
data_before_screen = data_all_conds.Raw;

%% =========================================================================
%% 5B. 2D GRID-SEARCH OPTIMIZATION (TARGETING FIRST 10 TRIALS PER CONDITION)
%% =========================================================================
cfg = struct('win', [0 1.0], 'base', [-1.0 -0.8], 'erp_win', [-0.2 0.8], 'm', 6);
iWin      = time_s >= cfg.win(1)     & time_s < cfg.win(2);
iBase     = time_s >= cfg.base(1)    & time_s < cfg.base(2);
iErpWin   = time_s >= cfg.erp_win(1) & time_s <= cfg.erp_win(2);
iZeroZone = time_s >= -0.200         & time_s <= 0.200; % Peristimulus baseline window [-0.2, 0.2] s
iPreStim  = time_s >= -0.200         & time_s <= 0.000; % Expectation ERP window [-0.2, 0.0] s

cond_keys     = {'P1', 'P2_500', 'P3_500', 'P2_2000', 'P3_missing'};
smooth_pts    = round(0.100 * 500); % 100 ms moving average for sustained drift
sat_thresh_uV = 100;                % Hard peak-amplitude ceiling: +/- 100 uV

% --- 1. PRE-COMPUTE TRIAL METRICS FOR ALL TARGET SUBJECTS & CONDITIONS ---
precomp = struct();
for t_idx = 1:length(target_subjs)
    s = target_subjs(t_idx);
    for c = 1:length(cond_keys)
        c_name = cond_keys{c};
        X_raw = data_all_conds.Raw.(c_name){s};
        [~, ~, nTr] = size(X_raw);
        
        X_bc_all       = zeros(size(X_raw));
        trial_max_amp  = zeros(nTr, 1);
        trial_c3_drift = zeros(nTr, 1);
        
        for tr = 1:nTr
            epoch = X_raw(:, :, tr);
            epoch_bc = epoch - mean(epoch(:, iBase), 2);
            X_bc_all(:, :, tr) = epoch_bc;
            
            % Bidirectional C3 Drift: Max |100ms-smoothed voltage| in [-0.2, 0.8] s
            c3_smooth = movmean(epoch_bc(chan_idx, iErpWin), smooth_pts);
            trial_c3_drift(tr) = max(abs(c3_smooth));
            
            % Peak Absolute Amplitude on C3 in [-0.2, 0.8] s
            trial_max_amp(tr) = max(abs(epoch_bc(chan_idx, iErpWin)));
        end
        
        precomp(s).(c_name).X_bc_all       = X_bc_all;
        precomp(s).(c_name).trial_c3_drift = trial_c3_drift;
        precomp(s).(c_name).trial_max_amp  = trial_max_amp;
        precomp(s).(c_name).c3_traces      = squeeze(X_bc_all(chan_idx, :, :)); % [nTime x nTr]
    end
end

% --- 2. 2D GRID SEARCH TARGETING FIRST 10 TRIALS OF EACH CONDITION ---
opt_mode    = 'per_subject';   % 'per_subject' or 'global'
drift_grid  = 16 : 1 : 50;     % Candidate drift thresholds: 16 to 50 uV
mad_grid    = 2.5 : 0.25 : 10; % Candidate MAD multipliers: 2.5 to 10.0
min_keep_p3 = 14;              % Minimum kept trials out of 30
ref_idx_10  = 1:10;            % First 10 trials as condition-specific expectation baseline

cost_surface   = inf(length(drift_grid), length(mad_grid), length(target_subjs));
null_diff_f10  = nan(length(drift_grid), length(mad_grid), length(target_subjs));
r500_diff_f10  = nan(length(drift_grid), length(mad_grid), length(target_subjs));
null_abs_val   = nan(length(drift_grid), length(mad_grid), length(target_subjs));
r500_abs_val   = nan(length(drift_grid), length(mad_grid), length(target_subjs));
rej_counts_mis = zeros(length(drift_grid), length(mad_grid), length(target_subjs));
rej_counts_500 = zeros(length(drift_grid), length(mad_grid), length(target_subjs));

f10_null_base = zeros(length(target_subjs), 1);
f10_r500_base = zeros(length(target_subjs), 1);

for t_idx = 1:length(target_subjs)
    s = target_subjs(t_idx);
    
    % Extract P3_missing (Rand Null) and P3_500 (Rand 500) metrics
    d_mis  = precomp(s).P3_missing.trial_c3_drift;
    a_mis  = precomp(s).P3_missing.trial_max_amp;
    tr_mis = precomp(s).P3_missing.c3_traces;
    
    d_500  = precomp(s).P3_500.trial_c3_drift;
    a_500  = precomp(s).P3_500.trial_max_amp;
    tr_500 = precomp(s).P3_500.c3_traces;
    
    % Compute condition-specific First-10-Trials Reference ERPs (excluding any >35 uV spike in 1:10)
    valid_10_mis = ref_idx_10(d_mis(ref_idx_10) <= 35 & a_mis(ref_idx_10) <= sat_thresh_uV);
    if isempty(valid_10_mis), valid_10_mis = ref_idx_10; end
    erp_mis_f10 = mean(tr_mis(:, valid_10_mis), 2, 'omitnan');
    
    valid_10_500 = ref_idx_10(d_500(ref_idx_10) <= 35 & a_500(ref_idx_10) <= sat_thresh_uV);
    if isempty(valid_10_500), valid_10_500 = ref_idx_10; end
    erp_500_f10 = mean(tr_500(:, valid_10_500), 2, 'omitnan');
    
    % Store First-10-Trials reference traces for plotting in Section 7
    precomp(s).P3_missing.erp_f10 = erp_mis_f10;
    precomp(s).P3_500.erp_f10     = erp_500_f10;
    
    f10_null_base(t_idx) = mean(erp_mis_f10(iPreStim));
    f10_r500_base(t_idx) = mean(erp_500_f10(iPreStim));
    
    for id = 1:length(drift_grid)
        d_lim = drift_grid(id);
        for im = 1:length(mad_grid)
            m_val = mad_grid(im);
            
            % --- Evaluate Rand Null (P3_missing) ---
            rob_d_mis = median(d_mis) + m_val * mad(d_mis, 1);
            rob_a_mis = median(a_mis) + m_val * mad(a_mis, 1);
            bad_mis   = (a_mis > sat_thresh_uV) | (a_mis > rob_a_mis) | ...
                        (d_mis > d_lim)         | (d_mis > rob_d_mis);
            
            % --- Evaluate Rand 500 (P3_500) ---
            rob_d_500 = median(d_500) + m_val * mad(d_500, 1);
            rob_a_500 = median(a_500) + m_val * mad(a_500, 1);
            bad_500   = (a_500 > sat_thresh_uV) | (a_500 > rob_a_500) | ...
                        (d_500 > d_lim)         | (d_500 > rob_d_500);
            
            if sum(~bad_500) < min_keep_p3
                [~, ord500] = sort(d_500, 'ascend');
                bad_500 = true(size(d_500));
                bad_500(ord500(1:min_keep_p3)) = false;
            end
            if sum(~bad_mis) < min_keep_p3
                continue; 
            end
            
            erp_mis_cand = mean(tr_mis(:, ~bad_mis), 2, 'omitnan');
            erp_500_cand = mean(tr_500(:, ~bad_500), 2, 'omitnan');
            
            % Objective: Match each condition's cleaned ERP to its OWN First-10-Trials expectation baseline
            % 1. Rand Null alignment to First 10 Trials ([-0.2, 0.0]s expectation + [-0.2, 0.2]s peristimulus + full trace)
            err_mis_pre  = abs(mean(erp_mis_cand(iPreStim))  - mean(erp_mis_f10(iPreStim)));
            err_mis_zero = abs(mean(erp_mis_cand(iZeroZone)) - mean(erp_mis_f10(iZeroZone)));
            err_mis_full = abs(mean(erp_mis_cand(iErpWin))   - mean(erp_mis_f10(iErpWin)));
            
            % 2. Rand 500 alignment to First 10 Trials ([-0.2, 0.0]s expectation + [-0.2, 0.2]s peristimulus)
            err_500_pre  = abs(mean(erp_500_cand(iPreStim))  - mean(erp_500_f10(iPreStim)));
            err_500_zero = abs(mean(erp_500_cand(iZeroZone)) - mean(erp_500_f10(iZeroZone)));
            
            % 3. Penalize keeping trials that exceed 30 uV so optimizer doesn't inflate drift_thresh_uV
            %    to cancel positive drifts with negative drifts
            max_kept_drift = max([d_mis(~bad_mis); d_500(~bad_500)]);
            leak_penalty   = 0.15 * max(0, max_kept_drift - 28);
            
            cost_surface(id, im, t_idx) = 1.2 * err_mis_pre + 1.0 * err_mis_zero + 0.5 * err_mis_full + ...
                                          1.0 * err_500_pre + 0.5 * err_500_zero + leak_penalty;
                                      
            null_diff_f10(id, im, t_idx)  = mean(erp_mis_cand(iPreStim)) - mean(erp_mis_f10(iPreStim));
            r500_diff_f10(id, im, t_idx)  = mean(erp_500_cand(iPreStim)) - mean(erp_500_f10(iPreStim));
            null_abs_val(id, im, t_idx)   = mean(erp_mis_cand(iPreStim));
            r500_abs_val(id, im, t_idx)   = mean(erp_500_cand(iPreStim));
            rej_counts_mis(id, im, t_idx) = sum(bad_mis);
            rej_counts_500(id, im, t_idx) = sum(bad_500);
        end
    end
end

% Determine optimal (drift_thresh_uV, mad_mult)
subj_opt_drift = zeros(length(target_subjs), 1);
subj_opt_mad   = zeros(length(target_subjs), 1);

joint_cost = sum(cost_surface, 3);
[~, lin_idx_glob] = min(joint_cost(:));
[id_glob, im_glob] = ind2sub(size(joint_cost), lin_idx_glob);
glob_drift_opt = drift_grid(id_glob);
glob_mad_opt   = mad_grid(im_glob);

fprintf('\n========================================================================================\n');
fprintf(' 2D GRID-SEARCH RESULTS (TARGETING FIRST 10 TRIALS BASELINE [-0.2, 0.0s] PER CONDITION)\n');
fprintf('========================================================================================\n');
fprintf('>>> GLOBAL OPTIMUM (Shared S3 & S14): drift_thresh_uV = %d uV | mad_mult = %.2f\n', ...
    glob_drift_opt, glob_mad_opt);

for t_idx = 1:length(target_subjs)
    s = target_subjs(t_idx);
    sub_cost = cost_surface(:, :, t_idx);
    [~, lin_idx_s] = min(sub_cost(:));
    [id_s, im_s] = ind2sub(size(sub_cost), lin_idx_s);
    
    if strcmpi(opt_mode, 'per_subject')
        subj_opt_drift(t_idx) = drift_grid(id_s);
        subj_opt_mad(t_idx)   = mad_grid(im_s);
        id_app = id_s; im_app = im_s;
    else
        subj_opt_drift(t_idx) = glob_drift_opt;
        subj_opt_mad(t_idx)   = glob_mad_opt;
        id_app = id_glob; im_app = im_glob;
    end
    
    fprintf('\n  * Subject %02d First-10-Trials Baseline [-0.2,0s] -> Rand Null: %+5.2f uV | Rand 500: %+5.2f uV\n', ...
        s, f10_null_base(t_idx), f10_r500_base(t_idx));
    fprintf('    - Applied (%s): drift_thresh_uV = %d uV | mad_mult = %.2f\n', ...
        upper(opt_mode), subj_opt_drift(t_idx), subj_opt_mad(t_idx));
    fprintf('    - Cleaned Rand Null [-0.2,0s]: %+5.2f uV (Diff vs First 10: %+5.2f uV | Rej: %d/30)\n', ...
        null_abs_val(id_app, im_app, t_idx), null_diff_f10(id_app, im_app, t_idx), rej_counts_mis(id_app, im_app, t_idx));
    fprintf('    - Cleaned Rand 500  [-0.2,0s]: %+5.2f uV (Diff vs First 10: %+5.2f uV | Rej: %d/30)\n', ...
        r500_abs_val(id_app, im_app, t_idx), r500_diff_f10(id_app, im_app, t_idx), rej_counts_500(id_app, im_app, t_idx));
end

% --- 3. APPLY OPTIMAL THRESHOLDS & EXPORT SCREENED DATASETS ---
screening_stats = struct();

for t_idx = 1:length(target_subjs)
    s = target_subjs(t_idx);
    drift_thresh_uV = subj_opt_drift(t_idx);
    mad_mult        = subj_opt_mad(t_idx);
    
    fprintf('\n--- SUBJECT %02d TRIAL SCREENING SUMMARY (drift_thresh = %d uV, mad_mult = %.2f) ---\n', ...
        s, drift_thresh_uV, mad_mult);
    fprintf('%-12s | %-10s | %-12s | %-12s | %-12s | %-10s\n', ...
        'Condition', 'Total Tr', 'Saturated', 'C3 Drift', 'Total Rej', 'Kept Tr');
    fprintf('--------------------------------------------------------------------------------\n');
    
    for c = 1:length(cond_keys)
        c_name = cond_keys{c};
        X_bc_all       = precomp(s).(c_name).X_bc_all;
        trial_c3_drift = precomp(s).(c_name).trial_c3_drift;
        trial_max_amp  = precomp(s).(c_name).trial_max_amp;
        nTr            = length(trial_c3_drift);
        
        rob_sat_lim   = median(trial_max_amp)  + mad_mult * mad(trial_max_amp, 1);
        rob_drift_lim = median(trial_c3_drift) + mad_mult * mad(trial_c3_drift, 1);
        
        is_saturated = (trial_max_amp > sat_thresh_uV)    | (trial_max_amp > rob_sat_lim);
        is_drifting  = (trial_c3_drift > drift_thresh_uV) | (trial_c3_drift > rob_drift_lim);
        bad_trials   = is_saturated | is_drifting;
        
        min_keep = min(min_keep_p3, nTr);
        if sum(~bad_trials) < min_keep
            artifact_score = (trial_c3_drift / (median(trial_c3_drift) + eps)) + ...
                             (trial_max_amp  / (median(trial_max_amp) + eps));
            [~, ord_clean] = sort(artifact_score, 'ascend');
            bad_trials = true(nTr, 1);
            bad_trials(ord_clean(1:min_keep)) = false;
        end
        
        X_clean = X_bc_all;
        X_clean(:, :, bad_trials) = NaN;
        data_all_conds.Raw.(c_name){s} = X_clean;
        
        screening_stats(s).(c_name).bad_trials      = bad_trials;
        screening_stats(s).(c_name).trial_max_amp   = trial_max_amp;
        screening_stats(s).(c_name).trial_c3_drift  = trial_c3_drift;
        screening_stats(s).(c_name).X_before        = X_bc_all;
        screening_stats(s).(c_name).X_after         = X_clean;
        screening_stats(s).(c_name).drift_thresh_uV = drift_thresh_uV;
        screening_stats(s).(c_name).mad_mult        = mad_mult;
        if isfield(precomp(s).(c_name), 'erp_f10')
            screening_stats(s).(c_name).erp_f10 = precomp(s).(c_name).erp_f10;
        end
        
        fprintf('%-12s | %-10d | %-12d | %-12d | %-12d | %-10d\n', ...
            c_name, nTr, sum(is_saturated), sum(is_drifting), sum(bad_trials), sum(~bad_trials));
    end
    
    if save_trial_cleaned_sets
        EEG_p1 = EEG_templates.P1{s};
        EEG_p1.data(:, :, 1:2:EEG_p1.trials) = single(data_all_conds.Raw.P1{s});
        out_p1 = fullfile(input_path, 'P1_cleaned');
        if ~exist(out_p1, 'dir'), mkdir(out_p1); end
        EEG_p1.icaact = []; % Clear mismatched activity matrix
        EEG_p1 = eeg_checkset(EEG_p1); % Recalculate internal consistency checks
        pop_saveset(EEG_p1, 'filename', EEG_p1.filename_orig, 'filepath', out_p1);
        
        EEG_p2 = EEG_templates.P2{s};
        EEG_p2.data(:, :, epoch_trials_p2_500)  = single(data_all_conds.Raw.P2_500{s});
        EEG_p2.data(:, :, epoch_trials_p2_2000) = single(data_all_conds.Raw.P2_2000{s});
        out_p2 = fullfile(input_path, 'P2_cleaned');
        if ~exist(out_p2, 'dir'), mkdir(out_p2); end
        EEG_p2.icaact = []; % Clear mismatched activity matrix
        EEG_p2 = eeg_checkset(EEG_p2); % Recalculate internal consistency checks
        pop_saveset(EEG_p2, 'filename', EEG_p2.filename_orig, 'filepath', out_p2);
        
        EEG_p3 = EEG_templates.P3{s};
        EEG_p3.data(:, :, epoch_trials_p3_500) = single(data_all_conds.Raw.P3_500{s});
        EEG_p3.data(:, :, epoch_trials_p3_mis) = single(data_all_conds.Raw.P3_missing{s});
        out_p3 = fullfile(input_path, 'P3_cleaned');
        if ~exist(out_p3, 'dir'), mkdir(out_p3); end
        EEG_p3.icaact = []; % Clear mismatched activity matrix
        EEG_p3 = eeg_checkset(EEG_p3); % Recalculate internal consistency checks
        pop_saveset(EEG_p3, 'filename', EEG_p3.filename_orig, 'filepath', out_p3);
    end
end

%% =========================================================================
%% 6. RE-CALCULATE DIAGNOSTIC METRICS AFTER TRIAL CLEANUP
%% =========================================================================
covCentered = @(X, b) (( (X - b) - mean(X - b, 2) ) * ( (X - b) - mean(X - b, 2) )') / size(X, 2);
V_ceil_all   = zeros(num_subjs, 1);
Base_Var_all = zeros(num_subjs, 5);

for s = valid_subjs
    raw_ref = data_all_conds.Raw.P1{s};
    half_n = floor(size(raw_ref, 3) / 2);
    erpA = mean(raw_ref(:, :, 1:half_n), 3, 'omitnan');
    erpB = mean(raw_ref(:, :, (half_n+1):(2*half_n)), 3, 'omitnan');
    
    bA = mean(erpA(:, iBase), 2, 'omitnan');
    bB = mean(erpB(:, iBase), 2, 'omitnan');
    C_A = covCentered(erpA(:, iWin), bA);
    [V, D] = eig(C_A, 'vector'); [~, ord] = sort(D, 'descend');
    Xi = V(:, ord(1:cfg.m));
    
    XB = erpB(:, iWin) - bB;
    V_ceil_all(s) = sum(sum((Xi' * XB).^2, 2) / nnz(iWin));
    
    for c = 1:5
        erp_all = mean(data_all_conds.Raw.(cond_keys{c}){s}, 3, 'omitnan');
        Base_Var_all(s, c) = mean(var(erp_all(:, iBase), 0, 2));
    end
end

normative_subjs = setdiff(valid_subjs, target_subjs);
median_Vceil = median(V_ceil_all(normative_subjs));

for s = target_subjs
    Vceil_ratio = V_ceil_all(s) / median_Vceil;
    base_drift_500 = Base_Var_all(s, 3) / (Base_Var_all(s, 1) + eps);
    base_drift_mis = Base_Var_all(s, 5) / (Base_Var_all(s, 1) + eps);
    
    fprintf('\n========================================================================\n');
    fprintf(' POST-CLEANUP DIAGNOSTIC RESULT FOR SUBJECT %d\n', s);
    fprintf('========================================================================\n');
    fprintf('1. Denominator (V_ceil) Ratio           : %6.2f%% of Cohort Median\n', Vceil_ratio * 100);
    fprintf('2. Baseline Noise Drift (P3_500 / P1)   : %6.2f x\n', base_drift_500);
    fprintf('3. Baseline Noise Drift (P3_mis / P1)   : %6.2f x\n', base_drift_mis);
end

%% =========================================================================
%% 7. PLOT TRIAL-LEVEL SCREENING & C3 ERP VERIFICATION (S3 & S14)
%% =========================================================================
for t_idx = 1:length(target_subjs)
    s = target_subjs(t_idx);
    curr_drift_thresh = screening_stats(s).P3_missing.drift_thresh_uV;
    curr_mad_mult     = screening_stats(s).P3_missing.mad_mult;
    
    fig_verify = figure('Position', [50, 50, 1700, 900], ...
        'Name', sprintf('Subject %d C3 Trial Screening & ERP Cleanup', s));
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    % Pre-calculate baseline-corrected ERPs for Panels 3 and 6
    erp_p1_clean = squeeze(mean(data_all_conds.Raw.P1{s}(chan_idx, :, :), 3, 'omitnan'));
    erp_500_bef  = squeeze(mean(screening_stats(s).P3_500.X_before(chan_idx, :, :), 3, 'omitnan'));
    erp_500_aft  = squeeze(mean(data_all_conds.Raw.P3_500{s}(chan_idx, :, :), 3, 'omitnan'));
    erp_500_f10  = screening_stats(s).P3_500.erp_f10;
    
    erp_mis_bef  = squeeze(mean(screening_stats(s).P3_missing.X_before(chan_idx, :, :), 3, 'omitnan'));
    erp_mis_aft  = squeeze(mean(data_all_conds.Raw.P3_missing{s}(chan_idx, :, :), 3, 'omitnan'));
    erp_mis_f10  = screening_stats(s).P3_missing.erp_f10;

    % --- PRE-CALCULATE SYNCHRONIZED Y-LIMITS FOR EACH COLUMN ---
    drift_v_500 = screening_stats(s).P3_500.trial_c3_drift;
    sat_v_500   = screening_stats(s).P3_500.trial_max_amp;
    drift_v_mis = screening_stats(s).P3_missing.trial_c3_drift;
    sat_v_mis   = screening_stats(s).P3_missing.trial_max_amp;
    col1_max    = max([drift_v_500(:); sat_v_500(:); drift_v_mis(:); sat_v_mis(:); sat_thresh_uV]);
    if isempty(col1_max) || isnan(col1_max) || col1_max <= 0, col1_max = 100; end
    ylim_col1   = [0, col1_max * 1.05];
    
    X_bef_500 = screening_stats(s).P3_500.X_before;
    X_bef_mis = screening_stats(s).P3_missing.X_before;
    tr_500_bc = squeeze(X_bef_500(chan_idx, :, :));
    tr_mis_bc = squeeze(X_bef_mis(chan_idx, :, :));
    
    col2_bound = max(abs([tr_500_bc(:); tr_mis_bc(:)]));
    if isempty(col2_bound) || isnan(col2_bound) || col2_bound <= 0, col2_bound = 100; end
    ylim_col2  = [-col2_bound * 1.05, col2_bound * 1.05];
    
    % Force all sliced erp arrays into strictly 1D column vectors before concatenating
    v_p1  = erp_p1_clean(iErpWin);
    v_5b  = erp_500_bef(iErpWin);
    v_5a  = erp_500_aft(iErpWin);
    v_5f  = erp_500_f10(iErpWin);
    v_mb  = erp_mis_bef(iErpWin);
    v_ma  = erp_mis_aft(iErpWin);
    v_mf  = erp_mis_f10(iErpWin);
    
    all_erps = [v_p1(:); v_5b(:); v_5a(:); v_5f(:); v_mb(:); v_ma(:); v_mf(:)];
    all_erps = all_erps(~isnan(all_erps));
    
    if isempty(all_erps)
        ylim_col3 = [-5, 5];
    else
        col3_min = min(all_erps);
        col3_max = max(all_erps);
        pad3     = (col3_max - col3_min) * 0.15;
        if pad3 == 0 || isnan(pad3), pad3 = 2; end
        ylim_col3 = [col3_min - pad3, col3_max + pad3];
    end
    if ylim_col3(1) >= ylim_col3(2)
        ylim_col3 = [ylim_col3(1) - 5, ylim_col3(2) + 5];
    end

    % ================= ROW 1: RAND 500 =================
    nexttile(1); hold on; set(gca, 'FontSize', 14);
    bad_tr_500 = screening_stats(s).P3_500.bad_trials;
    tr_idx_500 = 1:length(bad_tr_500);
    
    plot(tr_idx_500, drift_v_500, 'b-o', 'LineWidth', 1.5, 'DisplayName', '|C3 Drift| (\muV)');
    plot(tr_idx_500, sat_v_500, 'k-s', 'LineWidth', 1.5, 'DisplayName', '|C3 Peak Amp| (\muV)');
    if any(bad_tr_500)
        plot(tr_idx_500(bad_tr_500), sat_v_500(bad_tr_500), 'rx', 'MarkerSize', 10, 'LineWidth', 2.5, 'DisplayName', 'Rejected');
        plot(tr_idx_500(bad_tr_500), drift_v_500(bad_tr_500), 'rx', 'MarkerSize', 10, 'LineWidth', 2.5, 'HandleVisibility', 'off');
    end
    yline(curr_drift_thresh, 'b--', sprintf('Drift Thresh (%g \\muV)', curr_drift_thresh), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(sat_thresh_uV, 'r--', sprintf('Sat Thresh (%g \\muV)', sat_thresh_uV), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xlabel('Trial Index (Rand 500)'); ylabel('|Amplitude / Drift| (\muV)');
    title(sprintf('Rand 500 Screening (Rej: %d/%d)', sum(bad_tr_500), length(bad_tr_500)), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on; ylim(ylim_col1);
    
    nexttile(2); hold on; set(gca, 'FontSize', 14);
    for tr = 1:size(X_bef_500, 3)
        trace_c3 = tr_500_bc(:, tr);
        if bad_tr_500(tr)
            plot(time_s, trace_c3, 'Color', [0.9 0.2 0.2 0.6], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        else
            plot(time_s, trace_c3, 'Color', [0.6 0.6 0.6 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
        end
    end
    plot(nan, nan, 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'DisplayName', 'Kept');
    plot(nan, nan, 'Color', [0.9 0.2 0.2], 'LineWidth', 2, 'DisplayName', 'Rejected');
    yline(curr_drift_thresh, 'b--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(-curr_drift_thresh, 'b--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlim([-0.2, 0.8]); ylim(ylim_col2);
    xlabel('Time (s)'); ylabel('C3 Voltage (\muV)');
    title(sprintf('Rand 500 Traces (\\pm%g \\muV Bounds)', curr_drift_thresh), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on;
    
    nexttile(3); hold on; set(gca, 'FontSize', 14);
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(0, 'k:', 'HandleVisibility', 'off');
    plot(time_s, erp_500_bef, '--', 'Color', [0.9 0.6 0.1 0.4], 'LineWidth', 1.8, 'DisplayName', 'Rand 500 (Before)');
    plot(time_s, erp_p1_clean, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Cued P1 (Cleaned)');
    plot(time_s, erp_500_f10, ':', 'Color', [0.55 0.35 0.0], 'LineWidth', 2.2, 'DisplayName', 'Rand 500 (Trials 1-10 Ref)');
    plot(time_s, erp_500_aft, '-', 'Color', [0.9 0.6 0.1], 'LineWidth', 2.5, 'DisplayName', 'Rand 500 (Cleaned)');
    xlim([-0.2, 0.8]); ylim(ylim_col3);
    xlabel('Time (s)'); ylabel('Voltage (\muV)');
    title(sprintf('S%d C3 ERP: Rand 500', s), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on;

    % ================= ROW 2: RAND NULL =================
    nexttile(4); hold on; set(gca, 'FontSize', 14);
    bad_tr_mis = screening_stats(s).P3_missing.bad_trials;
    tr_idx_mis = 1:length(bad_tr_mis);
    
    plot(tr_idx_mis, drift_v_mis, 'b-o', 'LineWidth', 1.5, 'DisplayName', '|C3 Drift| (\muV)');
    plot(tr_idx_mis, sat_v_mis, 'k-s', 'LineWidth', 1.5, 'DisplayName', '|C3 Peak Amp| (\muV)');
    if any(bad_tr_mis)
        plot(tr_idx_mis(bad_tr_mis), sat_v_mis(bad_tr_mis), 'rx', 'MarkerSize', 10, 'LineWidth', 2.5, 'DisplayName', 'Rejected');
        plot(tr_idx_mis(bad_tr_mis), drift_v_mis(bad_tr_mis), 'rx', 'MarkerSize', 10, 'LineWidth', 2.5, 'HandleVisibility', 'off');
    end
    yline(curr_drift_thresh, 'b--', sprintf('Drift Thresh (%g \\muV)', curr_drift_thresh), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(sat_thresh_uV, 'r--', sprintf('Sat Thresh (%g \\muV)', sat_thresh_uV), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xlabel('Trial Index (Rand Null)'); ylabel('|Amplitude / Drift| (\muV)');
    title(sprintf('Rand Null Screening (Rej: %d/%d)', sum(bad_tr_mis), length(bad_tr_mis)), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on; ylim(ylim_col1);
    
    nexttile(5); hold on; set(gca, 'FontSize', 14);
    for tr = 1:size(X_bef_mis, 3)
        trace_c3 = tr_mis_bc(:, tr);
        if bad_tr_mis(tr)
            plot(time_s, trace_c3, 'Color', [0.9 0.2 0.2 0.6], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        else
            plot(time_s, trace_c3, 'Color', [0.6 0.6 0.6 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
        end
    end
    plot(nan, nan, 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'DisplayName', 'Kept');
    plot(nan, nan, 'Color', [0.9 0.2 0.2], 'LineWidth', 2, 'DisplayName', 'Rejected');
    yline(curr_drift_thresh, 'b--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(-curr_drift_thresh, 'b--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlim([-0.2, 0.8]); ylim(ylim_col2);
    xlabel('Time (s)'); ylabel('C3 Voltage (\muV)');
    title(sprintf('Rand Null Traces (\\pm%g \\muV Bounds)', curr_drift_thresh), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on;
    
    nexttile(6); hold on; set(gca, 'FontSize', 14);
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(0, 'k:', 'HandleVisibility', 'off');
    plot(time_s, erp_mis_bef, '--', 'Color', [0.2 0.7 0.9 0.4], 'LineWidth', 1.8, 'DisplayName', 'Rand Null (Before)');
    plot(time_s, erp_p1_clean, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Cued P1 (Cleaned)');
    plot(time_s, erp_mis_f10, ':', 'Color', [0.0 0.35 0.65], 'LineWidth', 2.2, 'DisplayName', 'Rand Null (Trials 1-10 Ref)');
    plot(time_s, erp_mis_aft, '-', 'Color', [0.2 0.7 0.9], 'LineWidth', 2.5, 'DisplayName', 'Rand Null (Cleaned)');
    xlim([-0.2, 0.8]); ylim(ylim_col3);
    xlabel('Time (s)'); ylabel('Voltage (\muV)');
    title(sprintf('S%d C3 ERP: Rand Null', s), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on;
    
    sgtitle(sprintf('Subject %d Trial Screening (Drift Thresh = %g \\muV, MAD Mult = %.2f)', ...
        s, curr_drift_thresh, curr_mad_mult), 'FontSize', 18, 'FontWeight', 'bold');
end