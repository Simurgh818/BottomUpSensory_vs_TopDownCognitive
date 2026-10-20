clear; clc; close all;

% --- Paths ---
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
    error('Unknown system: Cannot determine input and output paths.');
end

output_path = fullfile(base_output_path, 'ch-ch_cov');
if ~exist(output_path, 'dir')
    mkdir(output_path);
    fprintf('Created new output directory: %s\n', output_path);
end

screening_fig_dir = fullfile(output_path, 'Trial_Screening_Figures');
if ~exist(screening_fig_dir, 'dir')
    mkdir(screening_fig_dir);
end

conditions = {'BLA','BLT','P1','P2','P3'};
num_ch = 32;
fs = []; 
time_ms_eeg = [];
save_trial_cleaned_sets = true;

%% 0.1 Load and import file as EEGLAB matrix
data_all_conds = struct();
data_all_conds.Raw = struct();
data_all_conds.Alpha = struct();
data_all_conds.Beta = struct();
EEG_templates = struct();

% --- PRE-LOAD EXCEL INDICES ONCE ---
excel_p2 = fullfile(input_path, 'Indexes for P2.xlsx');
if exist(excel_p2, 'file')
    epoch_trials_p2_500  = readmatrix(excel_p2, 'Sheet', 'Audio onset with 500 ms tactile');
    epoch_trials_p2_2000 = readmatrix(excel_p2, 'Sheet', 'Audio onset with 2000 ms tactil');
    epoch_trials_p2_500  = epoch_trials_p2_500(~isnan(epoch_trials_p2_500));
    epoch_trials_p2_2000 = epoch_trials_p2_2000(~isnan(epoch_trials_p2_2000));
else
    error('Could not find: %s', excel_p2);
end

excel_p3 = fullfile(input_path, 'Indexes for P3.xlsx');
if exist(excel_p3, 'file')
    epoch_trials_p3_500 = readmatrix(excel_p3, 'Sheet', 'Audio onset with 500 ms tactile');
    epoch_trials_p3_mis = readmatrix(excel_p3, 'Sheet', 'Audio onset with missing tactil');
    epoch_trials_p3_500 = epoch_trials_p3_500(~isnan(epoch_trials_p3_500));
    epoch_trials_p3_mis = epoch_trials_p3_mis(~isnan(epoch_trials_p3_mis));
else
    error('Could not find: %s', excel_p3);
end

for c_idx = 1:length(conditions)
    cond_name = conditions{c_idx};
    in_dir = fullfile(input_path, cond_name);
    
    set_files = dir(fullfile(in_dir, '*.set'));
    names_sorted = sort(cellstr({set_files.name})');
    num_subjs_cond = length(names_sorted);
    
    if isempty(names_sorted)
        fprintf('No .set files for %s — skipping\n', cond_name);
        continue;
    end
    
    for s = 1:num_subjs_cond
        file_to_load = names_sorted{s};
        
        fprintf('Loading RAW Subj %d/%d (%s) for condition: %s\n', ...
            s, num_subjs_cond, file_to_load, cond_name);
        
        EEG = pop_loadset('filename', file_to_load, 'filepath', in_dir);
        
        if save_trial_cleaned_sets && ismember(cond_name, {'P1', 'P2', 'P3'})
            EEG_templates.(cond_name){s} = EEG;
            EEG_templates.(cond_name){s}.filename_orig = file_to_load;
        end
        
        if isempty(fs), fs = EEG.srate; end
        if isempty(time_ms_eeg)
            time_ms_eeg = linspace(0, 3.5, size(EEG.data, 2)); 
        end
        
        if ismember(cond_name, {'BLA', 'BLT', 'P1'})
            epoch_trials = 1:2:EEG.trials; % 60 trials for P1, BLA, BLT
            data_all_conds.Raw.(cond_name){s} = double(EEG.data(:, :, epoch_trials));
            
        elseif strcmp(cond_name, 'P2')
            data_all_conds.Raw.P2{s}      = double(EEG.data(:, :, 1:2:EEG.trials));
            data_all_conds.Raw.P2_500{s}  = double(EEG.data(:, :, epoch_trials_p2_500));
            data_all_conds.Raw.P2_2000{s} = double(EEG.data(:, :, epoch_trials_p2_2000));
            
        elseif strcmp(cond_name, 'P3')
            data_all_conds.Raw.P3{s}         = double(EEG.data(:, :, 1:2:EEG.trials));
            data_all_conds.Raw.P3_500{s}     = double(EEG.data(:, :, epoch_trials_p3_500));
            data_all_conds.Raw.P3_missing{s} = double(EEG.data(:, :, epoch_trials_p3_mis));
        end
    end
end

%% =========================================================================
%% SECTION 0.2: BANDPASS FILTERING (ALPHA, BETA, RATIO)
%% =========================================================================
disp('Applying Zero-Phase FIR Filters for Alpha and Beta...');

if isfield(EEG, 'srate') && ~isempty(EEG.srate) && isscalar(EEG.srate) && EEG.srate > 0
    fs = EEG.srate;
elseif exist('time_s', 'var') && length(time_s) > 1
    fs = round(1 / (time_s(2) - time_s(1)));
else
    fs = 500; 
end

fn = fs / 2;
fprintf('Filtering sample rate: %d Hz (Nyquist: %d Hz)\n', fs, fn);

ord_alpha = 3 * fix(fs / 8);
ord_beta  = 3 * fix(fs / 15);

b_alpha   = fir1(ord_alpha, [8 12] ./ fn, 'bandpass');
b_beta    = fir1(ord_beta,  [15 30] ./ fn, 'bandpass');

smooth_win = round(0.100 * fs);
smoothing_kernel = ones(smooth_win, 1) / smooth_win;

cond_names = fieldnames(data_all_conds.Raw);
for c = 1:length(cond_names)
    c_name = cond_names{c};
    for s = 1:length(data_all_conds.Raw.(c_name))
        raw_data = data_all_conds.Raw.(c_name){s};
        if isempty(raw_data), continue; end
        
        [nc, nt, ntr] = size(raw_data);
        
        X_perm = permute(raw_data, [2, 1, 3]);
        X_2D   = reshape(X_perm, nt, nc * ntr);
        
        X_alpha_2D = filtfilt(b_alpha, 1, double(X_2D));
        X_beta_2D  = filtfilt(b_beta, 1, double(X_2D));
        
        data_all_conds.Alpha.(c_name){s} = permute(reshape(X_alpha_2D, nt, nc, ntr), [2, 1, 3]);
        data_all_conds.Beta.(c_name){s}  = permute(reshape(X_beta_2D, nt, nc, ntr), [2, 1, 3]);
        
        env_alpha = abs(hilbert(X_alpha_2D));
        env_beta  = abs(hilbert(X_beta_2D));
        
        env_alpha_smooth = filtfilt(smoothing_kernel, 1, env_alpha);
        env_beta_smooth  = filtfilt(smoothing_kernel, 1, env_beta);
        
        ratio_2D = env_beta_smooth ./ (env_alpha_smooth + eps);
        data_all_conds.BetaAlphaRatio.(c_name){s} = permute(reshape(ratio_2D, nt, nc, ntr), [2, 1, 3]);
    end
end
disp('Filtering and Beta/Alpha Ratio computation complete.');

%% 1. Configuration Parameters
cfg = struct();
cfg.num_ch        = 32;
cfg.fs            = fs;                
cfg.tStart        = -1.000;            
cfg.win           = [0.000, 1.000];    
cfg.base          = [-1.000, -0.500];  % 500 ms pre-cue baseline [-1.0, -0.5 s]
cfg.nTrialMatch   = 20;                % Matched to 2/3 clean trial floor (20 trials)
cfg.nRep          = 50;                
cfg.m             = 6;                 
cfg.seed          = 20260910;
cfg.slide_win     = 0.100;
cfg.slide_step    = 0.025;

clean_name = @(c) strrep(strrep(strrep(strrep(strrep(strrep(strrep(c, ...
    'BLT', 'Tactile (BLT)'), 'BLA', 'Auditory (BLA)'), 'P1', 'Cued (P1)'), ...
    'P2_500', 'Unpred 500 (P2)'), 'P2_2000', 'Unpred 2000 (P2)'), ...
    'P3_500', 'Rand 500 (P3)'), 'P3_missing', 'Rand Null (P3)');

%% 2. Time Axis Alignment & Window Masking
nPnts = size(data_all_conds.Raw.P1{1}, 2); 
time_s = cfg.tStart + (0:(nPnts - 1)) / cfg.fs;
time_s = time_s(:)'; 

iWin   = time_s >= cfg.win(1)  & time_s < cfg.win(2);
iBase  = time_s >= cfg.base(1) & time_s < cfg.base(2);
T_win  = nnz(iWin);
t_eval = time_s(iWin);

t_starts    = time_s(1) : cfg.slide_step : (time_s(end) - cfg.slide_win);
t_centers   = t_starts + cfg.slide_win / 2;
num_windows = length(t_starts);

%% 3. Helper Functions (Unregularized Empirical Covariance)
covCentered = @(X, b) (( (X - b) - mean(X - b, 2) ) * ( (X - b) - mean(X - b, 2) )') / size(X, 2);

%% =========================================================================
%% SECTION 3B: UNIFIED FULL-EPOCH [-0.5, 2.4s] 32-CHANNEL SCREENING & 2/3 RULE
%% =========================================================================
num_all_subjs    = length(data_all_conds.Raw.P1);
all_target_subjs = 1:num_all_subjs;

chan_idx = 14; 
if exist('EEG', 'var') && isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs)
    c3_match = find(strcmpi({EEG.chanlocs.labels}, 'C3'), 1);
    if ~isempty(c3_match), chan_idx = c3_match; end
end
fprintf('\nExecuting Unified Full-Epoch [-0.5, 2.4s] 32-Channel Screening (95th Prctile, >=2/3 Keep Rule)...\n');

iErpWin   = time_s >= -0.200 & time_s <= 0.800; % Active ERP display window
iFullWin  = time_s >= -0.500 & time_s <= 2.400; % Full post-baseline trajectory window [-0.5, 2.4 s]
iZeroZone = time_s >= -0.200 & time_s <= 0.200;
iPreStim  = time_s >= -0.200 & time_s <= 0.000;

cond_keys_screen = {'P1', 'P2_500', 'P3_500', 'P2_2000', 'P3_missing', 'BLA', 'BLT'};
cond_keys_screen = cond_keys_screen(isfield(data_all_conds.Raw, cond_keys_screen));

smooth_pts    = round(0.100 * cfg.fs); % 100 ms moving average
sat_thresh_uV = 100;                   % 32-channel peak saturation ceiling (uV)

% --- UPDATED LENIENT GRID & 2/3 RETENTION FLOOR ---
drift_grid  = 25 : 1 : 60;       % Raised floor to 25 uV (up to 60 uV physical drift ceiling)
mad_grid    = 3.5 : 0.25 : 10.0; % Raised floor to 3.5 MAD
min_keep_30 = 20;                % Guarantee >= 20/30 (66.7%) kept trials
min_keep_60 = 40;                % Guarantee >= 40/60 (66.7%) kept trials

% --- 1. PRE-COMPUTE UNIFIED 32-CHANNEL FULL-EPOCH METRICS (95th Percentile) ---
precomp = struct();
for s = all_target_subjs
    for c = 1:length(cond_keys_screen)
        c_name = cond_keys_screen{c};
        X_raw = data_all_conds.Raw.(c_name){s};
        [~, ~, nTr] = size(X_raw);
        
        X_bc_all       = zeros(size(X_raw));
        trial_drift_32 = zeros(nTr, 1);
        trial_peak_32  = zeros(nTr, 1);
        
        for tr = 1:nTr
            epoch = X_raw(:, :, tr);
            epoch_bc = epoch - mean(epoch(:, iBase), 2);
            X_bc_all(:, :, tr) = epoch_bc;
            
            % 95th percentile across 32 channels in [-0.5, 2.4] s (requires >=2 channels with massive drift)
            ch_smooth_max = max(abs(movmean(epoch_bc(:, iFullWin), smooth_pts, 2)), [], 2);
            trial_drift_32(tr) = prctile(ch_smooth_max, 95);
            
            % 95th percentile peak amplitude across 32 channels in [-0.5, 2.4] s
            ch_peak_max = max(abs(epoch_bc(:, iFullWin)), [], 2);
            trial_peak_32(tr) = prctile(ch_peak_max, 95);
        end
        
        precomp(s).(c_name).X_bc_all       = X_bc_all;
        precomp(s).(c_name).trial_drift_32 = trial_drift_32;
        precomp(s).(c_name).trial_peak_32  = trial_peak_32;
        precomp(s).(c_name).c3_traces      = squeeze(X_bc_all(chan_idx, :, :));
    end
end

% --- 2. PER-SUBJECT, PER-CONDITION 2D GRID SEARCH (FIXED TRIAL CRITERION) ---
screening_stats = struct();
subj_excluded_2_3 = false(num_all_subjs, 1); % Tracks subjects that fail the 2/3 clean trial rule

tbl_Subj            = [];
tbl_Cond            = {};
tbl_TotalTr         = [];
tbl_KeptTr          = [];
tbl_RejTr           = [];
tbl_DriftThresh_32  = [];
tbl_MAD_Mult_32     = [];
tbl_MaxKept_Drift   = [];
tbl_MaxKept_MADR    = [];
tbl_Rej_By_Sat      = [];
tbl_Rej_By_Drift    = [];
tbl_MinKeep_Capped  = [];
tbl_Subj_Excluded   = [];
tbl_Outlier_Trials  = {};

fprintf('\n========================================================================================================================\n');
fprintf(' UNIFIED 32-CHANNEL FULL-EPOCH [-0.5, 2.4s] GRID-SEARCH RESULTS (95th Prctile, >=2/3 TRIAL RETENTION)\n');
fprintf('========================================================================================================================\n');

for s = all_target_subjs
    for c = 1:length(cond_keys_screen)
        c_name   = cond_keys_screen{c};
        X_bc_all = precomp(s).(c_name).X_bc_all;
        d_32     = precomp(s).(c_name).trial_drift_32;
        a_32     = precomp(s).(c_name).trial_peak_32;
        nTr      = length(d_32);
        
        if nTr >= 50
            min_keep = min_keep_60; % 40/60 (2/3)
        else
            min_keep = min(min_keep_30, ceil((2/3) * nTr)); % 20/30 (2/3)
        end
        
        best_cost  = inf;
        best_drift = 35;
        best_mad   = 4.0;
        
        for id = 1:length(drift_grid)
            d_lim = drift_grid(id);
            for im = 1:length(mad_grid)
                m_val = mad_grid(im);
                
                rob_d = median(d_32) + m_val * mad(d_32, 1);
                rob_a = median(a_32) + m_val * mad(a_32, 1);
                
                bad_cand = (a_32 > sat_thresh_uV) | (a_32 > rob_a) | ...
                           (d_32 > d_lim)         | (d_32 > rob_d);
                       
                if sum(~bad_cand) < min_keep
                    continue;
                end
                
                erp_32_cand = mean(X_bc_all(:, :, ~bad_cand), 3, 'omitnan');
                
                % Objective: Target 0 uV pre-stimulus baseline and peristimulus window across all trials
                err_pre  = mean(abs(mean(erp_32_cand(:, iPreStim), 2)));
                err_zero = mean(abs(mean(erp_32_cand(:, iZeroZone), 2)));
                err_full = mean(abs(mean(erp_32_cand(:, iFullWin), 2)));
                
                max_kept_d   = max(d_32(~bad_cand));
                leak_penalty = 0.15 * max(0, max_kept_d - 32) + 0.03 * sum(bad_cand);
                
                cost_val = 1.2 * err_pre + 1.0 * err_zero + 0.5 * err_full + leak_penalty;
                       
                if cost_val < best_cost
                    best_cost  = cost_val;
                    best_drift = d_lim;
                    best_mad   = m_val;
                end
            end
        end
        
        drift_thresh_uV = best_drift;
        mad_mult        = best_mad;
        
        mad_ratio_32  = (d_32 - median(d_32)) ./ (mad(d_32, 1) + eps);
        rob_sat_lim   = median(a_32) + mad_mult * mad(a_32, 1);
        rob_drift_lim = median(d_32) + mad_mult * mad(d_32, 1);
        
        is_saturated = (a_32 > sat_thresh_uV)    | (a_32 > rob_sat_lim);
        is_drifting  = (d_32 > drift_thresh_uV) | (d_32 > rob_drift_lim);
        bad_trials   = is_saturated | is_drifting;
        
        % Check if subject-condition fails to naturally retain >= 2/3 clean trials
        was_capped = false;
        if isinf(best_cost) || sum(~bad_trials) < min_keep
            was_capped = true;
            subj_excluded_2_3(s) = true; % Flag subject for exclusion (cannot keep 2/3 clean trials)
            
            artifact_score = (d_32 / (median(d_32) + eps)) + (a_32 / (median(a_32) + eps));
            [~, ord_clean] = sort(artifact_score, 'ascend');
            bad_trials = true(nTr, 1);
            bad_trials(ord_clean(1:min_keep)) = false;
        end
        
        X_clean = X_bc_all;
        X_clean(:, :, bad_trials) = NaN;
        data_all_conds.Raw.(c_name){s} = X_clean;
        
        if isfield(data_all_conds.Alpha, c_name)
            data_all_conds.Alpha.(c_name){s}(:, :, bad_trials) = NaN;
        end
        if isfield(data_all_conds.Beta, c_name)
            data_all_conds.Beta.(c_name){s}(:, :, bad_trials) = NaN;
        end
        if isfield(data_all_conds.BetaAlphaRatio, c_name)
            data_all_conds.BetaAlphaRatio.(c_name){s}(:, :, bad_trials) = NaN;
        end
        
        screening_stats(s).(c_name).bad_trials      = bad_trials;
        screening_stats(s).(c_name).trial_max_amp   = a_32;
        screening_stats(s).(c_name).trial_drift_32  = d_32;
        screening_stats(s).(c_name).mad_ratio_32    = mad_ratio_32;
        screening_stats(s).(c_name).X_before        = X_bc_all;
        screening_stats(s).(c_name).X_after         = X_clean;
        screening_stats(s).(c_name).drift_thresh_uV = drift_thresh_uV;
        screening_stats(s).(c_name).mad_mult        = mad_mult;
        
        rej_idx_list = find(bad_trials)';
        if isempty(rej_idx_list)
            rej_str = 'None';
        else
            rej_str = strjoin(arrayfun(@num2str, rej_idx_list, 'UniformOutput', false), ', ');
        end
        
        tbl_Subj(end+1, 1)           = s;
        tbl_Cond{end+1, 1}           = c_name;
        tbl_TotalTr(end+1, 1)        = nTr;
        tbl_KeptTr(end+1, 1)         = sum(~bad_trials);
        tbl_RejTr(end+1, 1)          = sum(bad_trials);
        tbl_DriftThresh_32(end+1, 1) = drift_thresh_uV;
        tbl_MAD_Mult_32(end+1, 1)    = mad_mult;
        tbl_MaxKept_Drift(end+1, 1)  = round(max([0; d_32(~bad_trials)]), 2);
        tbl_MaxKept_MADR(end+1, 1)   = round(max([0; mad_ratio_32(~bad_trials)]), 2);
        tbl_Rej_By_Sat(end+1, 1)     = sum(is_saturated & bad_trials);
        tbl_Rej_By_Drift(end+1, 1)   = sum(is_drifting & bad_trials);
        tbl_MinKeep_Capped(end+1, 1) = was_capped;
        tbl_Subj_Excluded(end+1, 1)  = false; % Populated after subject loop completes
        tbl_Outlier_Trials{end+1, 1} = rej_str;
    end
    
    % Mark all rows of this subject if they failed the 2/3 retention rule in any condition
    if subj_excluded_2_3(s)
        tbl_Subj_Excluded(tbl_Subj == s) = true;
        fprintf('  >>> WARNING: Subject %02d FAILED the >=2/3 clean trial retention rule and will be EXCLUDED from group analyses.\n', s);
    end
    
    if save_trial_cleaned_sets
        EEG_p1 = EEG_templates.P1{s};
        EEG_p1.data(:, :, 1:2:EEG_p1.trials) = single(data_all_conds.Raw.P1{s});
        out_p1 = fullfile(input_path, 'P1_cleaned');
        if ~exist(out_p1, 'dir'), mkdir(out_p1); end
        EEG_p1.icaact = []; EEG_p1 = eeg_checkset(EEG_p1);
        pop_saveset(EEG_p1, 'filename', EEG_p1.filename_orig, 'filepath', out_p1);
        
        EEG_p2 = EEG_templates.P2{s};
        EEG_p2.data(:, :, epoch_trials_p2_500)  = single(data_all_conds.Raw.P2_500{s});
        EEG_p2.data(:, :, epoch_trials_p2_2000) = single(data_all_conds.Raw.P2_2000{s});
        out_p2 = fullfile(input_path, 'P2_cleaned');
        if ~exist(out_p2, 'dir'), mkdir(out_p2); end
        EEG_p2.icaact = []; EEG_p2 = eeg_checkset(EEG_p2);
        pop_saveset(EEG_p2, 'filename', EEG_p2.filename_orig, 'filepath', out_p2);
        
        EEG_p3 = EEG_templates.P3{s};
        EEG_p3.data(:, :, epoch_trials_p3_500) = single(data_all_conds.Raw.P3_500{s});
        EEG_p3.data(:, :, epoch_trials_p3_mis) = single(data_all_conds.Raw.P3_missing{s});
        out_p3 = fullfile(input_path, 'P3_cleaned');
        if ~exist(out_p3, 'dir'), mkdir(out_p3); end
        EEG_p3.icaact = []; EEG_p3 = eeg_checkset(EEG_p3);
        pop_saveset(EEG_p3, 'filename', EEG_p3.filename_orig, 'filepath', out_p3);
    end
end

% --- 3. EXPORT OUTLIER TRIAL SUMMARY TABLE (.XLSX, .CSV, .TXT) ---
Outlier_Summary_Table = table(tbl_Subj, tbl_Cond, tbl_TotalTr, tbl_KeptTr, tbl_RejTr, ...
    tbl_DriftThresh_32, tbl_MAD_Mult_32, tbl_MaxKept_Drift, tbl_MaxKept_MADR, ...
    tbl_Rej_By_Sat, tbl_Rej_By_Drift, tbl_MinKeep_Capped, tbl_Subj_Excluded, tbl_Outlier_Trials, ...
    'VariableNames', {'Subject', 'Condition', 'Total_Trials', 'Kept_Trials', 'Rejected_Trials', ...
    'Drift_Thresh_32Ch_uV', 'MAD_Mult_32Ch', 'Max_Kept_32Ch_Drift_uV', 'Max_Kept_32Ch_MAD_Ratio', ...
    'Rej_By_Saturation', 'Rej_By_Drift', 'MinKeep_Floor_Triggered', 'Subject_Excluded', 'Outlier_Trial_Indices'});

csv_table_path  = fullfile(screening_fig_dir, 'Trial_Outlier_Summary_Table.csv');
xlsx_table_path = fullfile(screening_fig_dir, 'Trial_Outlier_Summary_Table.xlsx');
txt_table_path  = fullfile(screening_fig_dir, 'Trial_Outlier_Summary_Table.txt');

writetable(Outlier_Summary_Table, csv_table_path);
try writetable(Outlier_Summary_Table, xlsx_table_path); catch, end

fid_tbl = fopen(txt_table_path, 'w');
for out_dest = [1, fid_tbl]
    fprintf(out_dest, '\n==================================================================================================================================================================\n');
    fprintf(out_dest, ' UNIFIED 32-CHANNEL FULL-EPOCH [-0.5, 2.4s] TRIAL OUTLIER REPORT (95th Prctile, >=2/3 KEEP RULE)\n');
    fprintf(out_dest, '==================================================================================================================================================================\n');
    fprintf(out_dest, '%-4s | %-10s | %-7s | %-7s | %-18s | %-13s | %-13s | %-7s | %-8s | %s\n', ...
        'Subj', 'Condition', 'Kept/Tr', 'Rej (n)', '32Ch Thresh/MAD', 'MaxKeptDrift', 'MaxKeptMADR', 'Capped?', 'Excluded', 'Rejected Trial Indices');
    fprintf(out_dest, '------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
    for r = 1:height(Outlier_Summary_Table)
        cap_str = 'No'; if Outlier_Summary_Table.MinKeep_Floor_Triggered(r), cap_str = 'YES'; end
        exc_str = 'No'; if Outlier_Summary_Table.Subject_Excluded(r), exc_str = 'EXCLUDED'; end
        fprintf(out_dest, 'S%02d  | %-10s | %2d/%-4d | %2d      | %2d uV / %5.2f MAD  | %10.2f uV | %13.2f | %-7s | %-8s | %s\n', ...
            Outlier_Summary_Table.Subject(r), Outlier_Summary_Table.Condition{r}, ...
            Outlier_Summary_Table.Kept_Trials(r), Outlier_Summary_Table.Total_Trials(r), ...
            Outlier_Summary_Table.Rejected_Trials(r), ...
            Outlier_Summary_Table.Drift_Thresh_32Ch_uV(r), Outlier_Summary_Table.MAD_Mult_32Ch(r), ...
            Outlier_Summary_Table.Max_Kept_32Ch_Drift_uV(r), Outlier_Summary_Table.Max_Kept_32Ch_MAD_Ratio(r), ...
            cap_str, exc_str, Outlier_Summary_Table.Outlier_Trial_Indices{r});
    end
    fprintf(out_dest, '==================================================================================================================================================================\n\n');
end
fclose(fid_tbl);
fprintf('Saved Outlier Summary Tables to:\n  -> %s\n  -> %s\n  -> %s\n', csv_table_path, xlsx_table_path, txt_table_path);

% --- 4. SAVE 600 PPI PNG VERIFICATION FIGURES FOR ALL SUBJECTS ---
fprintf('\nGenerating and saving 600 ppi PNG Trial Screening figures for all %d subjects...\n', num_all_subjs);

for s = all_target_subjs
    d_thresh_500 = screening_stats(s).P3_500.drift_thresh_uV;
    mad_500      = screening_stats(s).P3_500.mad_mult;
    d_thresh_mis = screening_stats(s).P3_missing.drift_thresh_uV;
    mad_mis      = screening_stats(s).P3_missing.mad_mult;
    
    fig_verify = figure('Position', [50, 50, 1700, 900], 'Visible', 'off', ...
        'Name', sprintf('Subject %d 32-Ch Trial Screening & ERP Cleanup', s));
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    erp_p1_clean = squeeze(mean(data_all_conds.Raw.P1{s}(chan_idx, :, :), 3, 'omitnan'));
    erp_500_bef  = squeeze(mean(screening_stats(s).P3_500.X_before(chan_idx, :, :), 3, 'omitnan'));
    erp_500_aft  = squeeze(mean(data_all_conds.Raw.P3_500{s}(chan_idx, :, :), 3, 'omitnan'));
    
    erp_mis_bef  = squeeze(mean(screening_stats(s).P3_missing.X_before(chan_idx, :, :), 3, 'omitnan'));
    erp_mis_aft  = squeeze(mean(data_all_conds.Raw.P3_missing{s}(chan_idx, :, :), 3, 'omitnan'));

    drift_v_500 = screening_stats(s).P3_500.trial_drift_32;
    sat_v_500   = screening_stats(s).P3_500.trial_max_amp;
    drift_v_mis = screening_stats(s).P3_missing.trial_drift_32;
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
    
    v_p1 = erp_p1_clean(iErpWin);
    v_5b = erp_500_bef(iErpWin);
    v_5a = erp_500_aft(iErpWin);
    v_mb = erp_mis_bef(iErpWin);
    v_ma = erp_mis_aft(iErpWin);
    
    all_erps = [v_p1(:); v_5b(:); v_5a(:); v_mb(:); v_ma(:)];
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
    
    plot(tr_idx_500, drift_v_500, 'b-o', 'LineWidth', 1.5, 'DisplayName', '|32-Ch Drift| (\muV)');
    plot(tr_idx_500, sat_v_500, 'k-s', 'LineWidth', 1.5, 'DisplayName', '|32-Ch Peak| (\muV)');
    if any(bad_tr_500)
        plot(tr_idx_500(bad_tr_500), sat_v_500(bad_tr_500), 'rx', 'MarkerSize', 10, 'LineWidth', 2.5, 'DisplayName', 'Rejected');
        plot(tr_idx_500(bad_tr_500), drift_v_500(bad_tr_500), 'rx', 'MarkerSize', 10, 'LineWidth', 2.5, 'HandleVisibility', 'off');
    end
    yline(d_thresh_500, 'b--', sprintf('Drift Thresh (%g \\muV)', d_thresh_500), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(sat_thresh_uV, 'r--', sprintf('Sat Thresh (%g \\muV)', sat_thresh_uV), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xlabel('Trial Index (Rand 500)'); ylabel('32-Ch |Amplitude / Drift| (\muV)');
    title(sprintf('Rand 500 32-Ch Screening (Rej: %d/%d)', sum(bad_tr_500), length(bad_tr_500)), 'FontSize', 16, 'FontWeight', 'bold');
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
    yline(d_thresh_500, 'b--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(-d_thresh_500, 'b--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlim([-0.2, 0.8]); ylim(ylim_col2);
    xlabel('Time (s)'); ylabel('C3 Voltage (\muV)');
    title(sprintf('Rand 500 Traces (\\pm%g \\muV Bounds, %.2f MAD)', d_thresh_500, mad_500), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on;
    
    nexttile(3); hold on; set(gca, 'FontSize', 14);
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(0, 'k:', 'HandleVisibility', 'off');
    plot(time_s, erp_500_bef, '--', 'Color', [0.9 0.6 0.1 0.4], 'LineWidth', 1.8, 'DisplayName', 'Rand 500 (Before)');
    plot(time_s, erp_p1_clean, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Cued P1 (Cleaned)');
    plot(time_s, erp_500_aft, '-', 'Color', [0.9 0.6 0.1], 'LineWidth', 2.5, 'DisplayName', 'Rand 500 (Cleaned)');
    xlim([-0.2, 0.8]); ylim(ylim_col3);
    xlabel('Time (s)'); ylabel('Voltage (\muV)');
    title(sprintf('S%d C3 ERP: Rand 500', s), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on;

    % ================= ROW 2: RAND NULL =================
    nexttile(4); hold on; set(gca, 'FontSize', 14);
    bad_tr_mis = screening_stats(s).P3_missing.bad_trials;
    tr_idx_mis = 1:length(bad_tr_mis);
    
    plot(tr_idx_mis, drift_v_mis, 'b-o', 'LineWidth', 1.5, 'DisplayName', '|32-Ch Drift| (\muV)');
    plot(tr_idx_mis, sat_v_mis, 'k-s', 'LineWidth', 1.5, 'DisplayName', '|32-Ch Peak| (\muV)');
    if any(bad_tr_mis)
        plot(tr_idx_mis(bad_tr_mis), sat_v_mis(bad_tr_mis), 'rx', 'MarkerSize', 10, 'LineWidth', 2.5, 'DisplayName', 'Rejected');
        plot(tr_idx_mis(bad_tr_mis), drift_v_mis(bad_tr_mis), 'rx', 'MarkerSize', 10, 'LineWidth', 2.5, 'HandleVisibility', 'off');
    end
    yline(d_thresh_mis, 'b--', sprintf('Drift Thresh (%g \\muV)', d_thresh_mis), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(sat_thresh_uV, 'r--', sprintf('Sat Thresh (%g \\muV)', sat_thresh_uV), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xlabel('Trial Index (Rand Null)'); ylabel('32-Ch |Amplitude / Drift| (\muV)');
    title(sprintf('Rand Null 32-Ch Screening (Rej: %d/%d)', sum(bad_tr_mis), length(bad_tr_mis)), 'FontSize', 16, 'FontWeight', 'bold');
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
    yline(d_thresh_mis, 'b--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(-d_thresh_mis, 'b--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlim([-0.2, 0.8]); ylim(ylim_col2);
    xlabel('Time (s)'); ylabel('C3 Voltage (\muV)');
    title(sprintf('Rand Null Traces (\\pm%g \\muV Bounds, %.2f MAD)', d_thresh_mis, mad_mis), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on;
    
    nexttile(6); hold on; set(gca, 'FontSize', 14);
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(0, 'k:', 'HandleVisibility', 'off');
    plot(time_s, erp_mis_bef, '--', 'Color', [0.2 0.7 0.9 0.4], 'LineWidth', 1.8, 'DisplayName', 'Rand Null (Before)');
    plot(time_s, erp_p1_clean, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Cued P1 (Cleaned)');
    plot(time_s, erp_mis_aft, '-', 'Color', [0.2 0.7 0.9], 'LineWidth', 2.5, 'DisplayName', 'Rand Null (Cleaned)');
    xlim([-0.2, 0.8]); ylim(ylim_col3);
    xlabel('Time (s)'); ylabel('Voltage (\muV)');
    title(sprintf('S%d C3 ERP: Rand Null', s), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on;
    
    excl_tag = '';
    if subj_excluded_2_3(s), excl_tag = ' [EXCLUDED: <2/3 CLEAN TRIALS]'; end
    sgtitle(sprintf('Subject %d Unified 32-Ch Screening%s | Rand 500: %g \\muV, %.2f MAD | Rand Null: %g \\muV, %.2f MAD', ...
        s, excl_tag, d_thresh_500, mad_500, d_thresh_mis, mad_mis), 'FontSize', 18, 'FontWeight', 'bold');
    
    fig_png_path = fullfile(screening_fig_dir, sprintf('Subject_%02d_Trial_Screening_32Ch_600ppi.png', s));
    try
        exportgraphics(fig_verify, fig_png_path, 'Resolution', 600);
    catch
        print(fig_verify, fig_png_path, '-dpng', '-r600');
    end
    close(fig_verify);
    fprintf('  -> Saved 600 ppi PNG: %s\n', fig_png_path);
end
save(fullfile(screening_fig_dir, 'all_subjects_screening_stats.mat'), 'screening_stats', 'Outlier_Summary_Table', 'subj_excluded_2_3');

%% =========================================================================
%% MASTER LOOP: PROCESS FREQUENCY BANDS (EXCLUDING SUBJECTS FAILING 2/3 RULE)
%% =========================================================================
bands_to_process = {'Raw', 'Alpha', 'Beta', 'BetaAlphaRatio'};
cfg.ref_cond = 'P1';
cfg.group_A = {'P2_500', 'P3_500'};
cfg.group_B = {'P2_2000', 'P3_missing', 'BLA'};
cfg.all_conds = [{cfg.ref_cond}, cfg.group_A, cfg.group_B];

aggregate_spectra = struct();
aggregate_log_r   = struct();
aggregate_G       = struct();
aggregate_lat     = struct();
aggregate_sust    = struct();

subj_group_names = {'All_Subjects', 'Subjs_4_to_11', 'Other_Subjects'};
fig_tags         = {'All_Subjects', 'Cohort1', 'Cohort2'};

valid_retained_subjs = find(~subj_excluded_2_3(:))';
fprintf('\n======================================================\n');
fprintf(' RETAINED SUBJECTS FOR GROUP ANALYSIS (n=%d/%d): %s\n', ...
    length(valid_retained_subjs), num_all_subjs, mat2str(valid_retained_subjs));
if any(subj_excluded_2_3)
    fprintf(' EXCLUDED SUBJECTS (<2/3 Clean Trials): %s\n', mat2str(find(subj_excluded_2_3(:))'));
end
fprintf('======================================================\n');

for band_idx = 1:length(bands_to_process)
    band_name = bands_to_process{band_idx};
    current_data = data_all_conds.(band_name); 
    num_subjects = length(current_data.(cfg.ref_cond));
    
    fprintf('\n\n======================================================\n');
    fprintf('   BAND: %s | REFERENCE: %s (INDIVIDUAL COMPONENTS)\n', upper(band_name), cfg.ref_cond);
    fprintf('======================================================\n');
    
    group_dir = fullfile(output_path, 'group', band_name, cfg.ref_cond);
    if ~exist(group_dir, 'dir'), mkdir(group_dir); end

    Xi_all_subjs = nan(cfg.num_ch, cfg.m, num_subjects);
    lam_frac_all = nan(cfg.m, num_subjects);
    gap_all      = nan(cfg.m - 1, num_subjects);
    lam_full_all = nan(cfg.num_ch, num_subjects); 
    
    for s = valid_retained_subjs
        raw_ref = current_data.(cfg.ref_cond){s};
        valid_ref_idx = find(~isnan(squeeze(raw_ref(1, 1, :))));
        
        rng(cfg.seed + s);
        perm = valid_ref_idx(randperm(length(valid_ref_idx)));
        n_half = min(cfg.nTrialMatch, floor(length(valid_ref_idx) / 2));
        iA = perm(1:n_half);
        
        erpA = mean(raw_ref(:, :, iA), 3, 'omitnan');
        bA   = mean(erpA(:, iBase), 2, 'omitnan');
        C_A_cent = covCentered(erpA(:, iWin), bA);
        
        [V, D] = eig(C_A_cent, 'vector');
        [lam, ord] = sort(D, 'descend');
        V = V(:, ord);
        
        Xi_all_subjs(:, :, s) = V(:, 1:cfg.m);
        lam_frac_all(:, s)    = lam(1:cfg.m) / sum(lam);
        gap_all(:, s)         = (lam(1:cfg.m-1) - lam(2:cfg.m)) ./ lam(1:cfg.m-1);
        lam_full_all(:, s)    = lam / sum(lam); 
    end
    
    aggregate_spectra.(band_name) = mean(cumsum(lam_full_all(:, valid_retained_subjs), 1), 2, 'omitnan') * 100;
    
    Xi_aligned = Xi_all_subjs;
    ref_align_subj = valid_retained_subjs(1);
    for s = valid_retained_subjs(2:end)
        for k = 1:cfg.m
            if dot(Xi_aligned(:, k, s), Xi_aligned(:, k, ref_align_subj)) < 0
                Xi_aligned(:, k, s) = -Xi_aligned(:, k, s);
            end
        end
    end
    
    save(fullfile(group_dir, 'spectra.mat'), 'lam_frac_all', 'gap_all', 'valid_retained_subjs');
    
    figA = figure('Position', [100, 100, 1000, 400], 'Visible', 'off');
    subplot(1, 2, 1); plot(1:cfg.m, lam_frac_all(:, valid_retained_subjs), '-o', 'LineWidth', 1.5); grid on;
    xlabel('Direction Index'); ylabel('Fraction of Total Variance'); title(sprintf('Eigenvalue Spectra [%s|%s]', band_name, cfg.ref_cond));
    subplot(1, 2, 2); bar(mean(gap_all(:, valid_retained_subjs), 2, 'omitnan')); hold on;
    errorbar(1:cfg.m-1, mean(gap_all(:, valid_retained_subjs), 2, 'omitnan'), std(gap_all(:, valid_retained_subjs), 0, 2, 'omitnan'), 'k.', 'LineWidth', 1.2);
    grid on; xlabel('Gap Index'); ylabel('Relative Gap'); title('Mean Subspace Gaps');
    saveas(figA, fullfile(group_dir, 'Deliverable_A_Spectra.png')); close(figA);

    figC = figure('Position', [100, 100, 300 * cfg.m, 350], 'Visible', 'off');
    tiledlayout(1, cfg.m, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    n_ret = length(valid_retained_subjs);
    subj_subj_sim = nan(n_ret, n_ret, cfg.m);
    for k = 1:cfg.m
        for idx1 = 1:n_ret
            s1 = valid_retained_subjs(idx1);
            for idx2 = 1:n_ret
                s2 = valid_retained_subjs(idx2);
                subj_subj_sim(idx1, idx2, k) = abs(dot(Xi_aligned(:, k, s1), Xi_aligned(:, k, s2)));
            end
        end
        
        nexttile;
        imagesc(subj_subj_sim(:, :, k));
        colormap('jet'); clim([0 1]); axis square;
        title(sprintf('\\xi_%d', k), 'FontSize', 18, 'FontWeight', 'bold');
        
        if k == 1, ylabel('Retained Subject Index', 'FontSize', 12, 'FontWeight', 'bold'); end
        xlabel('Retained Subject Index', 'FontSize', 12, 'FontWeight', 'bold');
        xticks([1, n_ret]); yticks([1, n_ret]);
    end
    
    cb = colorbar; cb.Layout.Tile = 'east'; 
    cb.Label.String = 'Absolute Cosine Similarity'; cb.Label.FontSize = 14;
    sgtitle(sprintf('Cross-Subject Component Consistency (%s %s, n=%d)', band_name, cfg.ref_cond, n_ret), 'FontSize', 20, 'FontWeight', 'bold');
    
    save(fullfile(group_dir, 'subj_subj_similarity.mat'), 'subj_subj_sim', 'valid_retained_subjs');
    saveas(figC, fullfile(group_dir, 'Deliverable_C_SubjSubj_Sim.png')); close(figC);

    G_splits     = nan(length(cfg.all_conds), cfg.nRep, num_subjects);
    log_r_splits = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
    lat_splits   = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
    sust_splits  = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
    traj_splits  = nan(length(cfg.all_conds), cfg.m, num_windows, cfg.nRep, num_subjects);

    for subj = valid_retained_subjs
        subj_dir = fullfile(output_path, sprintf('subj%02d', subj), band_name, cfg.ref_cond);
        if ~exist(subj_dir, 'dir'), mkdir(subj_dir); end
        save(fullfile(subj_dir, 'cfg.mat'), 'cfg');
        
        Xi_stable = Xi_aligned(:, :, subj);
        
        for rep = 1:cfg.nRep
            rng(cfg.seed + rep + subj*1000); 
            
            raw_ref = current_data.(cfg.ref_cond){subj};
            valid_ref_idx = find(~isnan(squeeze(raw_ref(1, 1, :))));
            perm = valid_ref_idx(randperm(length(valid_ref_idx)));
            n_half = min(cfg.nTrialMatch, floor(length(valid_ref_idx) / 2));
            iA = perm(1:n_half);
            iB = perm(n_half + (1:n_half));
            
            erpA = mean(raw_ref(:, :, iA), 3, 'omitnan');
            bA   = mean(erpA(:, iBase), 2, 'omitnan');
            C_A_cent = covCentered(erpA(:, iWin), bA);
            
            [V, D] = eig(C_A_cent, 'vector');
            [lam, ord] = sort(D, 'descend');
            V = V(:, ord);
            Xi_rep = V(:, 1:cfg.m);
            lam_frac_rep = lam(1:cfg.m) / sum(lam);
            gap_rep = (lam(1:cfg.m-1) - lam(2:cfg.m)) ./ lam(1:cfg.m-1);
            
            for k = 1:cfg.m
                if dot(Xi_rep(:, k), Xi_stable(:, k)) < 0
                    Xi_rep(:, k) = -Xi_rep(:, k);
                end
            end
            
            basis_vars = struct('Xi', Xi_rep, 'lam', lam, 'lam_frac', lam_frac_rep, 'gap', gap_rep);
            
            erpB = mean(raw_ref(:, :, iB), 3, 'omitnan');
            bB   = mean(erpB(:, iBase), 2, 'omitnan');
            XB   = erpB(:, iWin) - bB;
            
            v_ceil = sum((Xi_rep' * XB).^2, 2) / T_win;
            V_ceil = sum(v_ceil);
            
            proj_B = Xi_rep' * XB; 
            [~, max_idx_B] = max(abs(proj_B), [], 2);
            t_lat_ceil = t_eval(max_idx_B); 
            
            split_vars = struct('iA', iA, 'iB', iB);
            proj_vars  = struct('v_ceil', v_ceil);
            
            for c = 1:length(cfg.all_conds)
                cond = cfg.all_conds{c};
                if ~isfield(current_data, cond), continue; end
                
                raw_cond = current_data.(cond){subj};
                if strcmp(cond, cfg.ref_cond)
                    Xk = XB; bk = bB; idx_k = iB;
                    erp_k = erpB; 
                else
                    valid_cond_idx = find(~isnan(squeeze(raw_cond(1, 1, :))));
                    n_pick = min(cfg.nTrialMatch, length(valid_cond_idx));
                    idx_k  = valid_cond_idx(randperm(length(valid_cond_idx), n_pick));
                    erp_k  = mean(raw_cond(:, :, idx_k), 3, 'omitnan');
                    bk     = mean(erp_k(:, iBase), 2, 'omitnan');
                    Xk     = erp_k(:, iWin) - bk;
                end
                
                split_vars.(cond) = idx_k;
                v_k = sum((Xi_rep' * Xk).^2, 2) / T_win;
                V_k = sum(v_k);
                G_k = V_k / V_ceil;
                r_k = v_k ./ v_ceil;
                log_r_k = log(r_k);
                
                proj_vars.(cond) = struct('v_k', v_k, 'G', G_k, 's_k', v_k/V_k, 'r', r_k, 'log_r', log_r_k);
                
                G_splits(c, rep, subj)        = G_k;
                log_r_splits(c, :, rep, subj) = log_r_k;
                
                proj_k = Xi_rep' * Xk; 
                [~, max_idx_k] = max(abs(proj_k), [], 2);
                lat_splits(c, :, rep, subj) = (t_eval(max_idx_k) - t_lat_ceil)' * 1000;
                
                xbar = mean(Xk, 2); Xd = Xk - xbar;
                v_sust = (Xi_rep' * xbar).^2;
                v_dyn  = sum((Xi_rep' * Xd).^2, 2) / T_win;
                sust_splits(c, :, rep, subj) = v_sust ./ (v_sust + v_dyn + eps);

                for w = 1:num_windows
                    w_mask = time_s >= t_starts(w) & time_s < (t_starts(w) + cfg.slide_win);
                    if nnz(w_mask) > 0
                        X_win = erp_k(:, w_mask) - bk;
                        v_win = sum((Xi_rep' * X_win).^2, 2) / nnz(w_mask);
                        traj_splits(c, :, w, rep, subj) = v_win;
                    end
                end
            end
            
            rep_str = sprintf('%02d', rep);
            save(fullfile(subj_dir, ['split_rep' rep_str '.mat']), '-struct', 'split_vars');
            save(fullfile(subj_dir, ['basis_rep' rep_str '.mat']), '-struct', 'basis_vars');
            save(fullfile(subj_dir, ['proj_rep' rep_str '.mat']), '-struct', 'proj_vars');
        end
    end

    % ------------------------------------------------------------------
    %% CHECKPOINT 1 (B): Cross-Subject Consistency Permutation Test
    % ------------------------------------------------------------------
    n_perms = 1000;
    num_pairs = (n_ret * (n_ret - 1)) / 2;
    true_mean_sim = zeros(cfg.m, 1);
    null_mean_sim = zeros(cfg.m, n_perms);
    p_vals        = zeros(cfg.m, 1);
    thresh_95     = zeros(cfg.m, 1);
    
    for k = 1:cfg.m
        pair_sims = zeros(num_pairs, 1);
        idx = 1;
        for i = 1:n_ret
            s1 = valid_retained_subjs(i);
            for j = (i+1):n_ret
                s2 = valid_retained_subjs(j);
                v1 = Xi_aligned(:, k, s1);
                v2 = Xi_aligned(:, k, s2);
                pair_sims(idx) = abs(dot(v1, v2)) / (norm(v1) * norm(v2));
                idx = idx + 1;
            end
        end
        true_mean_sim(k) = mean(pair_sims);
    end
    
    for p_idx = 1:n_perms
        for k = 1:cfg.m
            shuffled_Xi = zeros(cfg.num_ch, n_ret);
            for i = 1:n_ret
                s = valid_retained_subjs(i);
                shuff_idx = randperm(cfg.num_ch);
                shuffled_Xi(:, i) = Xi_aligned(shuff_idx, k, s);
            end
            
            pair_sims_null = zeros(num_pairs, 1);
            idx = 1;
            for i = 1:n_ret
                for j = (i+1):n_ret
                    v1 = shuffled_Xi(:, i);
                    v2 = shuffled_Xi(:, j);
                    pair_sims_null(idx) = abs(dot(v1, v2)) / (norm(v1) * norm(v2));
                    idx = idx + 1;
                end
            end
            null_mean_sim(k, p_idx) = mean(pair_sims_null);
        end
    end
    
    for k = 1:cfg.m
        p_vals(k) = (1 + sum(null_mean_sim(k, :) >= true_mean_sim(k))) / (1 + n_perms);
        thresh_95(k) = prctile(null_mean_sim(k, :), 95);
    end
    
    figC = figure('Position', [150, 150, 1400, 600], 'Name', sprintf('Permutation Test: %s', band_name), 'Visible', 'off');
    t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, sprintf('Cross-Subject Consistency: Null vs True (%s | %s, n=%d)', band_name, cfg.ref_cond, n_ret), 'FontSize', 18, 'FontWeight', 'bold');
    
    for k = 1:cfg.m
        nexttile; hold on;
        histogram(null_mean_sim(k, :), 30, 'Normalization', 'probability', 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
        xline(thresh_95(k), 'k--', 'LineWidth', 1.5, 'Label', '95% Thresh', 'LabelOrientation', 'horizontal', 'FontSize', 10);
        if p_vals(k) < 0.05, line_color = [0.8 0.1 0.1]; else, line_color = [0.2 0.2 0.2]; end
        xline(true_mean_sim(k), '-', 'Color', line_color, 'LineWidth', 2.5, 'Label', 'True Sim', 'LabelOrientation', 'horizontal', 'FontSize', 10);
        title(sprintf('\\xi_%d (p = %.3f)', k, p_vals(k)), 'FontSize', 16, 'FontWeight', 'bold');
        xlabel('Mean Abs. Cosine Similarity', 'FontSize', 12);
        if ismember(k, [1, 4]), ylabel('Probability', 'FontSize', 12); end
        xlim([0, max(max(null_mean_sim(k, :)), true_mean_sim(k)) * 1.1]);
        grid on; box on;
    end
    save_file_perm = fullfile(group_dir, 'Deliverable_C_Permutation_Test.png');
    try saveas(figC, save_file_perm); catch, end
    close(figC);

    % ------------------------------------------------------------------
    %% CHECKPOINT 2 & 3: COHORT DEFINITIONS (EXCLUDING FAILED SUBJECTS)
    % ------------------------------------------------------------------
    subj_mean_log_r = squeeze(mean(log_r_splits, 3, 'omitnan')); 
    subj_mean_traj  = squeeze(mean(traj_splits, 4, 'omitnan'));  
    
    all_subjs = valid_retained_subjs;                     % Retained All Subjects (e.g., n=13 without S14)
    g1_subjs  = intersect(4:11, valid_retained_subjs);    % Cohort A: Subjs 4 to 11 (n=8)
    g2_subjs  = setdiff(valid_retained_subjs, g1_subjs);  % Cohort B: Other Retained Subjects (e.g., n=5 without S14)
    
    subj_groups = {all_subjs, g1_subjs, g2_subjs};
    num_groups  = length(subj_groups);
    
    time_mask_zoom = (t_centers >= -0.1) & (t_centers <= 1.0);
    time_mask_full = (t_centers >= -0.15); 
    
    global_ylims_full = zeros(cfg.m, 2); 
    global_ylims_zoom = zeros(cfg.m, 2); 
    global_clims      = zeros(cfg.m, 1);
    all_log_r         = [];
    
    for k = 1:cfg.m
        all_up_full = []; all_dn_full = [];
        all_up_zoom = []; all_dn_zoom = [];
        max_c = 0; 
        
        for g_idx = 1:num_groups
            grp_s = subj_groups{g_idx};
            if isempty(grp_s), continue; end
            
            Xi_grp = mean(Xi_aligned(:, :, grp_s), 3, 'omitnan');
            Xi_grp = Xi_grp ./ vecnorm(Xi_grp);
            
            for c_idx = 1:length(cfg.all_conds)
                cond = cfg.all_conds{c_idx};
                g_val = squeeze(subj_mean_traj(c_idx, k, :, grp_s));
                if length(grp_s) == 1
                    g_val = g_val(:); g_m = g_val; g_s = zeros(size(g_val));
                else
                    g_m = mean(g_val, 2, 'omitnan');
                    g_s = std(g_val, 0, 2, 'omitnan') ./ sqrt(length(grp_s));
                end
                
                all_up_full = [all_up_full; g_m(time_mask_full) + g_s(time_mask_full)];
                all_dn_full = [all_dn_full; g_m(time_mask_full) - g_s(time_mask_full)];
                all_up_zoom = [all_up_zoom; g_m(time_mask_zoom) + g_s(time_mask_zoom)];
                all_dn_zoom = [all_dn_zoom; g_m(time_mask_zoom) - g_s(time_mask_zoom)];
                
                erp_grp = mean(cat(3, current_data.(cond){grp_s}), 3, 'omitnan');
                max_c = max(max_c, max(abs(Xi_grp(:, k)' * erp_grp)) * 0.5);
                
                if k == 1
                    g_lr_m = mean(subj_mean_log_r(c_idx, :, grp_s), 3, 'omitnan');
                    g_lr_s = std(subj_mean_log_r(c_idx, :, grp_s), 0, 3, 'omitnan') ./ sqrt(length(grp_s));
                    all_log_r = [all_log_r; g_lr_m(:) + g_lr_s(:); g_lr_m(:) - g_lr_s(:)];
                end
            end
        end
        
        if isempty(all_up_full), all_up_full = 1; all_dn_full = 0; end
        y_up_f = max(all_up_full); y_dn_f = min(all_dn_full);
        pad_f = (y_up_f - y_dn_f) * 0.05; if pad_f == 0 || isnan(pad_f), pad_f = 1; end
        global_ylims_full(k, 1) = max(0, y_dn_f - pad_f);
        global_ylims_full(k, 2) = y_up_f + pad_f;
        
        if isempty(all_up_zoom), all_up_zoom = 1; all_dn_zoom = 0; end
        y_up_z = max(all_up_zoom); y_dn_z = min(all_dn_zoom);
        pad_z = (y_up_z - y_dn_z) * 0.05; if pad_z == 0 || isnan(pad_z), pad_z = 1; end
        global_ylims_zoom(k, 1) = max(0, y_dn_z - pad_z);
        global_ylims_zoom(k, 2) = y_up_z + pad_z;
        
        if max_c == 0, max_c = 1; end 
        global_clims(k) = max_c;
    end
    
    if isempty(all_log_r), all_log_r = [-1 1]; end
    cp2_y_min = min(-1.5, floor(min(all_log_r) * 1.15));
    cp2_y_max = max( 1.5, ceil(max(all_log_r) * 1.15));

    for g = 1:num_groups
        curr_subjs = subj_groups{g};
        curr_n = length(curr_subjs);
        g_name = subj_group_names{g};
        if curr_n == 0, continue; end
        
        grp_mean_log_r = mean(subj_mean_log_r(:, :, curr_subjs), 3, 'omitnan');      
        grp_se_log_r   = std(subj_mean_log_r(:, :, curr_subjs), 0, 3, 'omitnan') ./ sqrt(curr_n);
        
        aggregate_log_r.(band_name).cohort(g).mean = grp_mean_log_r;
        aggregate_log_r.(band_name).cohort(g).se   = grp_se_log_r;
        aggregate_log_r.(band_name).cohort(g).subj_data = subj_mean_log_r(:, :, curr_subjs);

        subj_mean_G = squeeze(mean(G_splits(:, :, curr_subjs), 2, 'omitnan')); 
        if curr_n == 1, subj_mean_G = subj_mean_G(:); end
        aggregate_G.(band_name).cohort(g).mean = mean(subj_mean_G, 2, 'omitnan');
        aggregate_G.(band_name).cohort(g).se   = std(subj_mean_G, 0, 2, 'omitnan') ./ sqrt(curr_n);
        aggregate_G.(band_name).cohort(g).subj_data = subj_mean_G;
        
        subj_mean_lat = squeeze(mean(mean(lat_splits(:, :, :, curr_subjs), 3, 'omitnan'), 2, 'omitnan')); 
        if curr_n == 1, subj_mean_lat = subj_mean_lat(:); end
        lat_dir_mean = squeeze(mean(lat_splits(:, :, :, curr_subjs), 3, 'omitnan'));
        if curr_n == 1, lat_dir_mean = reshape(lat_dir_mean, [length(cfg.all_conds), cfg.m, 1]); end
        
        aggregate_lat.(band_name).cohort(g).mean = mean(lat_dir_mean, 3, 'omitnan');
        aggregate_lat.(band_name).cohort(g).se   = std(lat_dir_mean, 0, 3, 'omitnan') ./ sqrt(curr_n);
        aggregate_lat.(band_name).cohort(g).subj_data = subj_mean_lat; 
        
        subj_mean_sust = squeeze(mean(mean(sust_splits(:, :, :, curr_subjs), 3, 'omitnan'), 2, 'omitnan')) * 100;
        if curr_n == 1, subj_mean_sust = subj_mean_sust(:); end
        aggregate_sust.(band_name).cohort(g).mean = mean(subj_mean_sust, 2, 'omitnan');
        aggregate_sust.(band_name).cohort(g).se   = std(subj_mean_sust, 0, 2, 'omitnan') ./ sqrt(curr_n);
        aggregate_sust.(band_name).cohort(g).subj_data = subj_mean_sust; 
        
        % CHECKPOINT 3a TRAJECTORIES (COHORT-FITTED Y-LIMITS WITH 20% HEADROOM)
        grp_mean_traj = squeeze(mean(subj_mean_traj(:, :, :, curr_subjs), 4, 'omitnan')); 
        grp_sem_traj  = squeeze(std(subj_mean_traj(:, :, :, curr_subjs), 0, 4, 'omitnan')) ./ sqrt(curr_n);
        
        sig_components = find(p_vals < 0.05)';
        if isempty(sig_components), sig_components = 1; end
        num_sig = length(sig_components);
        
        figChk3 = figure('Position', [100, 100, 1600, max(600, num_sig * 500)], 'Visible', 'off');
        tiledlayout(num_sig, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        c_idx_ref = find(strcmp(cfg.all_conds, cfg.ref_cond));
        
        for k_idx = 1:num_sig
            k = sig_components(k_idx);
            row_mean_k = grp_mean_traj(:, k, time_mask_full);
            row_sem_k  = grp_sem_traj(:, k, time_mask_full);
            row_upper  = row_mean_k + row_sem_k;
            row_lower  = row_mean_k - row_sem_k;
            
            cp3a_min_y = max(0, min(row_lower(:)) * 0.90);
            cp3a_max_y = max(row_upper(:)) * 1.20;
            if isempty(cp3a_max_y) || isnan(cp3a_max_y) || cp3a_max_y <= cp3a_min_y
                cp3a_min_y = 0; cp3a_max_y = 10;
            end
            
            nexttile; hold on; set(gca, 'FontSize', 22);
            patch([cfg.win(1) cfg.win(2) cfg.win(2) cfg.win(1)], [cp3a_min_y cp3a_min_y cp3a_max_y cp3a_max_y], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
            xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
            
            if ~isempty(c_idx_ref)
                y_val_ref = reshape(grp_mean_traj(c_idx_ref, k, :), 1, []); y_sem_ref = reshape(grp_sem_traj(c_idx_ref, k, :), 1, []);
                fill([t_centers, fliplr(t_centers)], max(0, [y_val_ref + y_sem_ref, fliplr(y_val_ref - y_sem_ref)]), [0.5 0.5 0.5], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(t_centers, y_val_ref, 'k-', 'LineWidth', 2, 'DisplayName', sprintf('Reference (%s)', clean_name(cfg.ref_cond)));
            end
            
            colors_A = {[0.850 0.325 0.098], [0.929 0.694 0.125], [0.494 0.184 0.556]};
            for i = 1:length(cfg.group_A)
                c_idx = find(strcmp(cfg.all_conds, cfg.group_A{i}));
                if isempty(c_idx), continue; end
                y_val = reshape(grp_mean_traj(c_idx, k, :), 1, []); y_sem = reshape(grp_sem_traj(c_idx, k, :), 1, []);
                col_idx = mod(i-1, 3) + 1;
                fill([t_centers, fliplr(t_centers)], max(0, [y_val + y_sem, fliplr(y_val - y_sem)]), colors_A{col_idx}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(t_centers, y_val, '-', 'Color', colors_A{col_idx}, 'LineWidth', 2, 'DisplayName', clean_name(cfg.group_A{i}));
            end
            grid on; xlabel('Time (s)', 'FontSize', 20); ylabel('Projected Variance', 'FontSize', 20);
            title(sprintf('Group A: \\xi_%d Trajectory', k), 'FontSize', 24); 
            if k_idx == 1, lgdA = legend('Location', 'best'); lgdA.FontSize = 18; end
            xlim([t_centers(1), t_centers(end)]); ylim([cp3a_min_y, cp3a_max_y]);
            
            nexttile; hold on; set(gca, 'FontSize', 22);
            patch([cfg.win(1) cfg.win(2) cfg.win(2) cfg.win(1)], [cp3a_min_y cp3a_min_y cp3a_max_y cp3a_max_y], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
            xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
            
            if ~isempty(c_idx_ref)
                fill([t_centers, fliplr(t_centers)], max(0, [y_val_ref + y_sem_ref, fliplr(y_val_ref - y_sem_ref)]), [0.5 0.5 0.5], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(t_centers, y_val_ref, 'k-', 'LineWidth', 2, 'DisplayName', sprintf('Reference (%s)', clean_name(cfg.ref_cond)));
            end
            
            colors_B = {[0.466 0.674 0.188], [0.301 0.745 0.933], [0.494 0.184 0.556]};
            for i = 1:length(cfg.group_B)
                c_idx = find(strcmp(cfg.all_conds, cfg.group_B{i}));
                if isempty(c_idx), continue; end
                y_val = reshape(grp_mean_traj(c_idx, k, :), 1, []); y_sem = reshape(grp_sem_traj(c_idx, k, :), 1, []);
                col_idx = mod(i-1, length(colors_B)) + 1;
                fill([t_centers, fliplr(t_centers)], max(0, [y_val + y_sem, fliplr(y_val - y_sem)]), colors_B{col_idx}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(t_centers, y_val, '-', 'Color', colors_B{col_idx}, 'LineWidth', 2, 'DisplayName', clean_name(cfg.group_B{i}));
            end
            grid on; xlabel('Time (s)', 'FontSize', 20); ylabel('Projected Variance', 'FontSize', 20);
            title(sprintf('Group B: \\xi_%d Trajectory', k), 'FontSize', 24); 
            if k_idx == 1, lgdB = legend('Location', 'best'); lgdB.FontSize = 18; end
            xlim([t_centers(1), t_centers(end)]); ylim([cp3a_min_y, cp3a_max_y]);
        end
        sgtitle(sprintf('Subspace Trajectories: %s (%s, n=%d)', band_name, strrep(g_name, '_', ' '), curr_n), 'FontSize', 24, 'FontWeight', 'bold');
        
        save_base_chk3 = fullfile(group_dir, sprintf('Checkpoint3a_Trajectories_%s', g_name));
        try saveas(figChk3, [save_base_chk3, '.svg']); saveas(figChk3, [save_base_chk3, '.png']); catch, end
        close(figChk3);
    end 
end

disp('All analyses completed successfully with [-1.0 to -0.5 s] baseline, 95th percentile 32-Ch guard, and >=2/3 trial retention rule.');