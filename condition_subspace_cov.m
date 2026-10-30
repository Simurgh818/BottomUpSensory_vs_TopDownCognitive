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
    fs = 512; 
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

%% =========================================================================
%% SECTION 1: SYSTEM CONFIGURATION (SET TO 5 STATISTICALLY VALIDATED DIRECTIONS)
%% =========================================================================
cfg = struct();
cfg.num_ch        = 32;
cfg.fs            = fs;                
cfg.tStart        = -1.000;            
cfg.win           = [0.000, 1.000];    
cfg.base          = [-1.000, -0.500];  % 500 ms pre-cue baseline window [-1.0, -0.5 s]
cfg.nTrialMatch   = 20;                % Matched clean trial subsampling floor
cfg.nRep          = 50;                % Split-half Monte Carlo resampling repetitions
cfg.m             = 5;                 % TRUNCATED TO 5 VALIDATED MANIFOLD DIRECTIONS (p < 0.05)
cfg.seed          = 20260910;
cfg.slide_win     = 0.100;             % 100 ms sliding trajectory window
cfg.slide_step    = 0.025;             % 25 ms step size
cfg.basis_mode    = 'universal';       % Concatenated Universal Basis

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
%% SECTION 3B: FULL-EPOCH [-0.5, 2.4s] SCREENING (DRIFT + PEAK + PASER 6ms DERIVATIVE)
%% =========================================================================
num_all_subjs    = length(data_all_conds.Raw.P1);
all_target_subjs = 1:num_all_subjs;

% Locate C3, C4, F3, T7, and T8 indices from EEG.chanlocs
chan_c3 = 14; chan_c4 = 18; chan_f3 = 4; chan_t7 = 7; chan_t8 = 24;
non_ocular_chans = 1:cfg.num_ch;
if exist('EEG', 'var') && isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs)
    all_labels = {EEG.chanlocs.labels};
    m_c3 = find(strcmpi(all_labels, 'C3'), 1);
    m_c4 = find(strcmpi(all_labels, 'C4'), 1);
    m_f3 = find(strcmpi(all_labels, 'F3'), 1);
    m_t7 = find(ismember(upper(all_labels), {'T7', 'T3'}), 1);
    m_t8 = find(ismember(upper(all_labels), {'T8', 'T4'}), 1);
    if ~isempty(m_c3), chan_c3 = m_c3; end
    if ~isempty(m_c4), chan_c4 = m_c4; end
    if ~isempty(m_f3), chan_f3 = m_f3; end
    if ~isempty(m_t7), chan_t7 = m_t7; end
    if ~isempty(m_t8), chan_t8 = m_t8; end
    
    ocular_mask = ismember(upper(all_labels), {'FP1', 'FP2', 'FPZ', 'AF3', 'AF4', 'F7', 'F8'});
    if any(~ocular_mask), non_ocular_chans = find(~ocular_mask); end
end
chan_idx  = chan_c3; 
core_roi  = unique([chan_c3, chan_c4, chan_f3, chan_t7, chan_t8], 'stable');

fprintf('\nExecuting Screening on {{C3, C4, F3, T7, T8} + Non-Ocular 85th Prctile} (35:50 uV + PASER 6ms Derivative)...\n');

iErpWin   = time_s >= -0.200 & time_s <= 0.800; 
iCritWin  = time_s >= -0.500 & time_s <= 1.000; % Active task window [-0.5, 1.0 s]
iFullWin  = time_s >= -0.500 & time_s <= 2.400; % Full epoch [-0.5, 2.4 s]
iZeroZone = time_s >= -0.200 & time_s <= 0.200;
iPreStim  = time_s >= -0.200 & time_s <= 0.000;

cond_keys_screen = {'P1', 'P2_500', 'P3_500', 'P2_2000', 'P3_missing', 'BLA', 'BLT'};
cond_keys_screen = cond_keys_screen(isfield(data_all_conds.Raw, cond_keys_screen));

smooth_pts      = round(0.100 * cfg.fs);        % 100 ms moving-average window (50 samples)
deriv_pts       = max(1, round(0.006 * cfg.fs));% PASER Step 2.2: 6 ms derivative step (3 samples at 500 Hz)
sat_thresh_uV   = 100;                          % Peak saturation ceiling (uV)
deriv_thresh_uV = 60;                           % Physical 6 ms fast-step ceiling (60 uV / 6 ms)
paser_mad_floor = 6.0;                          % PASER Step 2.2: 6 MAD upper threshold for fast derivative
ext_roi         = core_roi;

% --- CONSERVATIVE GRID SEARCH & >=80% RETENTION PARAMETERS ---
drift_grid  = 35 : 1 : 50;       % Updated drift range: 35 to 50 uV
mad_grid    = 4.5 : 0.25 : 10.0; % Conservative MAD floor: 4.5 to 10.0 MAD
min_keep_30 = 24;                % Guarantee >= 24/30 (80.0%) kept trials (excludes Subj 3 & Subj 14)
min_keep_60 = 48;                % Guarantee >= 48/60 (80.0%) kept trials (excludes Subj 3 & Subj 14)

% --- 1. PRE-COMPUTE TRIAL METRICS (DRIFT, PEAK, AND PASER 6ms DERIVATIVE) ---
precomp = struct();
for s = all_target_subjs
    for c = 1:length(cond_keys_screen)
        c_name = cond_keys_screen{c};
        X_raw = data_all_conds.Raw.(c_name){s};
        [~, ~, nTr] = size(X_raw);
        
        X_bc_all       = zeros(size(X_raw));
        trial_drift_32 = zeros(nTr, 1);
        trial_peak_32  = zeros(nTr, 1);
        trial_deriv_32 = zeros(nTr, 1);
        
        for tr = 1:nTr
            epoch = X_raw(:, :, tr);
            epoch_bc = epoch - mean(epoch(:, iBase), 2);
            X_bc_all(:, :, tr) = epoch_bc;
            
            % 1. Smoothed slow drift (100 ms movmean)
            ch_smooth_crit = max(abs(movmean(epoch_bc(:, iCritWin), smooth_pts, 2)), [], 2);
            crit_drift     = max(max(ch_smooth_crit(core_roi)), prctile(ch_smooth_crit(non_ocular_chans), 85));
            
            ch_smooth_full = max(abs(movmean(epoch_bc(:, iFullWin), smooth_pts, 2)), [], 2);
            late_med_drift = 0.75 * max(max(ch_smooth_full(core_roi)), prctile(ch_smooth_full(non_ocular_chans), 85));
            trial_drift_32(tr) = max(crit_drift, late_med_drift);
            
            % 2. Peak amplitude (100 uV active / 200 uV full-epoch C.R.A.P.)
            ch_peak_crit = max(abs(epoch_bc(:, iCritWin)), [], 2);
            ch_peak_full = max(abs(epoch_bc(:, iFullWin)), [], 2);
            crit_peak    = max(max(ch_peak_crit(core_roi)), prctile(ch_peak_crit(non_ocular_chans), 85));
            crap_peak    = 0.5 * max(ch_peak_full(ext_roi)); 
            trial_peak_32(tr) = max(crit_peak, crap_peak);
            
            % 3. PASER Step 2.2: Fast 6 ms derivative |V(t + 6ms) - V(t)|
            E_crit        = epoch_bc(:, iCritWin);
            ch_deriv_crit = max(abs(E_crit(:, (1+deriv_pts):end) - E_crit(:, 1:(end-deriv_pts))), [], 2);
            trial_deriv_32(tr) = max(max(ch_deriv_crit(core_roi)), prctile(ch_deriv_crit(non_ocular_chans), 85));
        end
        
        precomp(s).(c_name).X_bc_all       = X_bc_all;
        precomp(s).(c_name).trial_drift_32 = trial_drift_32;
        precomp(s).(c_name).trial_peak_32  = trial_peak_32;
        precomp(s).(c_name).trial_deriv_32 = trial_deriv_32;
        precomp(s).(c_name).c3_traces      = squeeze(X_bc_all(chan_idx, :, :));
    end
end

% --- 2. PER-SUBJECT, PER-CONDITION 2D GRID SEARCH ---
screening_stats = struct();
subj_excluded_2_3 = false(num_all_subjs, 1);

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
tbl_Rej_By_Deriv    = [];
tbl_MinKeep_Capped  = [];
tbl_Subj_Excluded   = [];
tbl_Outlier_Trials  = {};

for s = all_target_subjs
    for c = 1:length(cond_keys_screen)
        c_name   = cond_keys_screen{c};
        X_bc_all = precomp(s).(c_name).X_bc_all;
        d_32     = precomp(s).(c_name).trial_drift_32;
        a_32     = precomp(s).(c_name).trial_peak_32;
        v_32     = precomp(s).(c_name).trial_deriv_32;
        nTr      = length(d_32);
        
        if nTr >= 50
            min_keep = min_keep_60; 
        else
            min_keep = min_keep_30; 
        end
        
        best_cost  = inf;
        best_drift = max(drift_grid);
        best_mad   = min(mad_grid);
        
        for id = 1:length(drift_grid)
            d_lim = drift_grid(id);
            for im = 1:length(mad_grid)
                m_val = mad_grid(im);
                
                rob_d = median(d_32) + m_val * mad(d_32, 1);
                rob_a = median(a_32) + m_val * mad(a_32, 1);
                rob_v = median(v_32) + max(paser_mad_floor, m_val) * mad(v_32, 1);
                
                bad_cand = (a_32 > sat_thresh_uV)   | (a_32 > rob_a) | ...
                           (d_32 > d_lim)           | (d_32 > rob_d) | ...
                           (v_32 > deriv_thresh_uV) | (v_32 > rob_v);
                       
                if sum(~bad_cand) < min_keep
                    continue;
                end
                
                erp_32_cand = mean(X_bc_all(:, :, ~bad_cand), 3, 'omitnan');
                
                err_pre  = mean(abs(mean(erp_32_cand(core_roi, iPreStim), 2)));
                err_zero = mean(abs(mean(erp_32_cand(core_roi, iZeroZone), 2)));
                err_full = median(abs(mean(erp_32_cand(:, iCritWin), 2)));
                rej_penalty = 0.15 * sum(bad_cand);
                
                cost_val = 1.2 * err_pre + 1.0 * err_zero + 0.5 * err_full + rej_penalty;
                       
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
        rob_deriv_lim = median(v_32) + max(paser_mad_floor, mad_mult) * mad(v_32, 1);
        
        is_saturated = (a_32 > sat_thresh_uV)    | (a_32 > rob_sat_lim);
        is_drifting  = (d_32 > drift_thresh_uV)  | (d_32 > rob_drift_lim);
        is_fast_step = (v_32 > deriv_thresh_uV)  | (v_32 > rob_deriv_lim);
        bad_trials   = is_saturated | is_drifting | is_fast_step;
        
        was_capped = false;
        if isinf(best_cost) || sum(~bad_trials) < min_keep
            was_capped = true;
            if ismember(c_name, {'P1', 'P2_500', 'P3_500', 'P2_2000', 'P3_missing', 'BLA'})
                subj_excluded_2_3(s) = true;
            end
            
            artifact_score = (d_32 / (median(d_32) + eps)) + ...
                             (a_32 / (median(a_32) + eps)) + ...
                             (v_32 / (median(v_32) + eps));
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
        screening_stats(s).(c_name).trial_deriv_32  = v_32;
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
        tbl_Rej_By_Deriv(end+1, 1)   = sum(is_fast_step & bad_trials);
        tbl_MinKeep_Capped(end+1, 1) = was_capped;
        tbl_Subj_Excluded(end+1, 1)  = false;
        tbl_Outlier_Trials{end+1, 1} = rej_str;
    end
    
    if subj_excluded_2_3(s)
        tbl_Subj_Excluded(tbl_Subj == s) = true;
        fprintf('  >>> WARNING: Subject %02d FAILED the >=80%% clean trial retention rule and will be EXCLUDED from group analyses.\n', s);
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

% --- 3. EXPORT OUTLIER TRIAL SUMMARY TABLE ---
Outlier_Summary_Table = table(tbl_Subj, tbl_Cond, tbl_TotalTr, tbl_KeptTr, tbl_RejTr, ...
    tbl_DriftThresh_32, tbl_MAD_Mult_32, tbl_MaxKept_Drift, tbl_MaxKept_MADR, ...
    tbl_Rej_By_Sat, tbl_Rej_By_Drift, tbl_Rej_By_Deriv, tbl_MinKeep_Capped, tbl_Subj_Excluded, tbl_Outlier_Trials, ...
    'VariableNames', {'Subject', 'Condition', 'Total_Trials', 'Kept_Trials', 'Rejected_Trials', ...
    'Drift_Thresh_32Ch_uV', 'MAD_Mult_32Ch', 'Max_Kept_32Ch_Drift_uV', 'Max_Kept_32Ch_MAD_Ratio', ...
    'Rej_By_Saturation', 'Rej_By_Drift', 'Rej_By_Derivative', 'MinKeep_Floor_Triggered', 'Subject_Excluded', 'Outlier_Trial_Indices'});

csv_table_path  = fullfile(screening_fig_dir, 'Trial_Outlier_Summary_Table.csv');
xlsx_table_path = fullfile(screening_fig_dir, 'Trial_Outlier_Summary_Table.xlsx');
txt_table_path  = fullfile(screening_fig_dir, 'Trial_Outlier_Summary_Table.txt');

writetable(Outlier_Summary_Table, csv_table_path);
try writetable(Outlier_Summary_Table, xlsx_table_path); catch, end

fid_tbl = fopen(txt_table_path, 'w');
for out_dest = [1, fid_tbl]
    fprintf(out_dest, '\n============================================================================================================================================================================\n');
    fprintf(out_dest, ' FULL-EPOCH [-0.5, 2.4s] TRIAL OUTLIER REPORT (35:50 uV + PASER 6ms DERIVATIVE, >=80%% KEEP RULE)\n');
    fprintf(out_dest, '============================================================================================================================================================================\n');
    fprintf(out_dest, '%-4s | %-10s | %-7s | %-7s | %-18s | %-13s | %-13s | %-9s | %-7s | %-8s | %s\n', ...
        'Subj', 'Condition', 'Kept/Tr', 'Rej (n)', 'Thresh / MAD', 'MaxKeptDrift', 'MaxKeptMADR', 'RejDeriv', 'Capped?', 'Excluded', 'Rejected Trial Indices');
    fprintf(out_dest, '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
    for r = 1:height(Outlier_Summary_Table)
        cap_str = 'No'; if Outlier_Summary_Table.MinKeep_Floor_Triggered(r), cap_str = 'YES'; end
        exc_str = 'No'; if Outlier_Summary_Table.Subject_Excluded(r), exc_str = 'EXCLUDED'; end
        fprintf(out_dest, 'S%02d  | %-10s | %2d/%-4d | %2d      | %2d uV / %5.2f MAD  | %10.2f uV | %13.2f | %9d | %-7s | %-8s | %s\n', ...
            Outlier_Summary_Table.Subject(r), Outlier_Summary_Table.Condition{r}, ...
            Outlier_Summary_Table.Kept_Trials(r), Outlier_Summary_Table.Total_Trials(r), ...
            Outlier_Summary_Table.Rejected_Trials(r), ...
            Outlier_Summary_Table.Drift_Thresh_32Ch_uV(r), Outlier_Summary_Table.MAD_Mult_32Ch(r), ...
            Outlier_Summary_Table.Max_Kept_32Ch_Drift_uV(r), Outlier_Summary_Table.Max_Kept_32Ch_MAD_Ratio(r), ...
            Outlier_Summary_Table.Rej_By_Derivative(r), ...
            cap_str, exc_str, Outlier_Summary_Table.Outlier_Trial_Indices{r});
    end
    fprintf(out_dest, '============================================================================================================================================================================\n\n');
end
fclose(fid_tbl);

% --- 4. SAVE 600 PPI PNG VERIFICATION FIGURES (WITH TOP TITLE INSIDE TILEDLAYOUT) ---
fprintf('\nGenerating and saving 600 ppi PNG Trial Screening figures for all %d subjects...\n', num_all_subjs);

for s = all_target_subjs
    d_thresh_500 = screening_stats(s).P3_500.drift_thresh_uV;
    mad_500      = screening_stats(s).P3_500.mad_mult;
    d_thresh_mis = screening_stats(s).P3_missing.drift_thresh_uV;
    mad_mis      = screening_stats(s).P3_missing.mad_mult;
    
    fig_verify = figure('Position', [50, 50, 1700, 950], 'Visible', 'off', 'Color', 'w', ...
        'Name', sprintf('Subject %d Trial Screening & ERP Cleanup', s));
    
    t_verify = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    erp_p1_clean = squeeze(mean(data_all_conds.Raw.P1{s}(chan_idx, :, :), 3, 'omitnan'));
    erp_500_bef  = squeeze(mean(screening_stats(s).P3_500.X_before(chan_idx, :, :), 3, 'omitnan'));
    erp_500_aft  = squeeze(mean(data_all_conds.Raw.P3_500{s}(chan_idx, :, :), 3, 'omitnan'));
    
    erp_mis_bef  = squeeze(mean(screening_stats(s).P3_missing.X_before(chan_idx, :, :), 3, 'omitnan'));
    erp_mis_aft  = squeeze(mean(data_all_conds.Raw.P3_missing{s}(chan_idx, :, :), 3, 'omitnan'));

    drift_v_500 = screening_stats(s).P3_500.trial_drift_32;
    sat_v_500   = screening_stats(s).P3_500.trial_max_amp;
    deriv_v_500 = screening_stats(s).P3_500.trial_deriv_32;
    drift_v_mis = screening_stats(s).P3_missing.trial_drift_32;
    sat_v_mis   = screening_stats(s).P3_missing.trial_max_amp;
    deriv_v_mis = screening_stats(s).P3_missing.trial_deriv_32;
    col1_max    = max([drift_v_500(:); sat_v_500(:); deriv_v_500(:); drift_v_mis(:); sat_v_mis(:); deriv_v_mis(:); sat_thresh_uV]);
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
    nexttile(t_verify, 1); hold on; set(gca, 'FontSize', 14);
    bad_tr_500 = screening_stats(s).P3_500.bad_trials;
    tr_idx_500 = 1:length(bad_tr_500);
    
    plot(tr_idx_500, drift_v_500, 'b-o', 'LineWidth', 1.5, 'DisplayName', '|ROI+Scalp Drift| (\muV)');
    plot(tr_idx_500, sat_v_500, 'k-s', 'LineWidth', 1.5, 'DisplayName', '|ROI+Scalp Peak| (\muV)');
    plot(tr_idx_500, deriv_v_500, 'g-^', 'Color', [0.15 0.65 0.25], 'LineWidth', 1.4, 'DisplayName', '|6ms Step Deriv| (\muV)');
    if any(bad_tr_500)
        plot(tr_idx_500(bad_tr_500), sat_v_500(bad_tr_500), 'rx', 'MarkerSize', 10, 'LineWidth', 2.5, 'DisplayName', 'Rejected');
        plot(tr_idx_500(bad_tr_500), drift_v_500(bad_tr_500), 'rx', 'MarkerSize', 10, 'LineWidth', 2.5, 'HandleVisibility', 'off');
    end
    yline(d_thresh_500, 'b--', sprintf('Drift Thresh (%g \\muV)', d_thresh_500), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(sat_thresh_uV, 'r--', sprintf('Sat Thresh (%g \\muV)', sat_thresh_uV), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xlabel('Trial Index (Rand 500)'); ylabel('|Amplitude / Drift| (\muV)');
    title(sprintf('Rand 500 Screening (Rej: %d/%d)', sum(bad_tr_500), length(bad_tr_500)), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on; ylim(ylim_col1);
    
    nexttile(t_verify, 2); hold on; set(gca, 'FontSize', 14);
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
    yline(0, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(d_thresh_500, 'b--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(-d_thresh_500, 'b--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlim([-0.5, 0.8]); ylim(ylim_col2);
    xlabel('Time (s)'); ylabel('C3 Voltage (\muV)');
    title(sprintf('Rand 500 Traces (\\pm%g \\muV Bounds, %.2f MAD)', d_thresh_500, mad_500), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on;
    
    nexttile(t_verify, 3); hold on; set(gca, 'FontSize', 14);
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
    nexttile(t_verify, 4); hold on; set(gca, 'FontSize', 14);
    bad_tr_mis = screening_stats(s).P3_missing.bad_trials;
    tr_idx_mis = 1:length(bad_tr_mis);
    
    plot(tr_idx_mis, drift_v_mis, 'b-o', 'LineWidth', 1.5, 'DisplayName', '|ROI+Scalp Drift| (\muV)');
    plot(tr_idx_mis, sat_v_mis, 'k-s', 'LineWidth', 1.5, 'DisplayName', '|ROI+Scalp Peak| (\muV)');
    plot(tr_idx_mis, deriv_v_mis, 'g-^', 'Color', [0.15 0.65 0.25], 'LineWidth', 1.4, 'DisplayName', '|6ms Step Deriv| (\muV)');
    if any(bad_tr_mis)
        plot(tr_idx_mis(bad_tr_mis), sat_v_mis(bad_tr_mis), 'rx', 'MarkerSize', 10, 'LineWidth', 2.5, 'DisplayName', 'Rejected');
        plot(tr_idx_mis(bad_tr_mis), drift_v_mis(bad_tr_mis), 'rx', 'MarkerSize', 10, 'LineWidth', 2.5, 'HandleVisibility', 'off');
    end
    yline(d_thresh_mis, 'b--', sprintf('Drift Thresh (%g \\muV)', d_thresh_mis), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(sat_thresh_uV, 'r--', sprintf('Sat Thresh (%g \\muV)', sat_thresh_uV), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xlabel('Trial Index (Rand Null)'); ylabel('|Amplitude / Drift| (\muV)');
    title(sprintf('Rand Null Screening (Rej: %d/%d)', sum(bad_tr_mis), length(bad_tr_mis)), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on; ylim(ylim_col1);
    
    nexttile(t_verify, 5); hold on; set(gca, 'FontSize', 14);
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
    yline(0, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(d_thresh_mis, 'b--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(-d_thresh_mis, 'b--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlim([-0.5, 0.8]); ylim(ylim_col2);
    xlabel('Time (s)'); ylabel('C3 Voltage (\muV)');
    title(sprintf('Rand Null Traces (\\pm%g \\muV Bounds, %.2f MAD)', d_thresh_mis, mad_mis), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on;
    
    nexttile(t_verify, 6); hold on; set(gca, 'FontSize', 14);
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(0, 'k:', 'HandleVisibility', 'off');
    plot(time_s, erp_mis_bef, '--', 'Color', [0.2 0.7 0.9 0.4], 'LineWidth', 1.8, 'DisplayName', 'Rand Null (Before)');
    plot(time_s, erp_p1_clean, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Cued P1 (Cleaned)');
    plot(time_s, erp_mis_aft, '-', 'Color', [0.2 0.7 0.9], 'LineWidth', 2.5, 'DisplayName', 'Rand Null (Cleaned)');
    xlim([-0.2, 0.8]); ylim(ylim_col3);
    xlabel('Time (s)'); ylabel('Voltage (\muV)');
    title(sprintf('S%d C3 ERP: Rand Null', s), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on;
    
    % Attach top title directly to t_verify so exportgraphics NEVER crops it out
    excl_tag = '';
    if subj_excluded_2_3(s), excl_tag = ' [EXCLUDED: <80% CLEAN TRIALS]'; end
    title(t_verify, sprintf('Subject %d Trial Screening%s (Rand 500: %g \\muV, %.2f MAD | Rand Null: %g \\muV, %.2f MAD)', ...
        s, excl_tag, d_thresh_500, mad_500, d_thresh_mis, mad_mis), 'FontSize', 20, 'FontWeight', 'bold');
    
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
%% MASTER LOOP: 5-DIMENSIONAL CONCATENATED UNIVERSAL BASIS & UNBIASED PROJECTION
%% =========================================================================
bands_to_process = {'Raw', 'Alpha', 'Beta', 'BetaAlphaRatio'};
cfg.ref_cond = 'P1';
cfg.group_A = {'P2_500', 'P3_500'};
cfg.group_B = {'P2_2000', 'P3_missing'}; % BLA excluded from Group B
cfg.all_conds = [{cfg.ref_cond}, cfg.group_A, cfg.group_B];

aggregate_spectra = struct();
aggregate_log_r   = struct();
aggregate_G       = struct();
aggregate_lat     = struct();
aggregate_sust    = struct();

% Single pooled cohort (All Retained Subjects, excluding S3 and S14)
subj_group_names = {'All_Subjects'};
fig_tags         = {'All_Subjects'};

subj_excluded_2_3(intersect([3, 14], 1:num_all_subjs)) = true;
valid_retained_subjs = find(~subj_excluded_2_3(:))';

fprintf('\n======================================================\n');
fprintf(' RETAINED SUBJECTS FOR GROUP ANALYSIS (n=%d/%d): %s\n', ...
    length(valid_retained_subjs), num_all_subjs, mat2str(valid_retained_subjs));
if any(subj_excluded_2_3)
    fprintf(' EXCLUDED SUBJECTS (<80%% Clean Trials): %s\n', mat2str(find(subj_excluded_2_3(:))'));
end
fprintf('======================================================\n');

% Close stale file handles and initialize safe text logging
fclose('all');
summary_log_file = fullfile(output_path, 'Summary_Tables_All_Bands.txt');
fid_sum = fopen(summary_log_file, 'w');
if fid_sum < 0
    warning('Could not open %s. Printing to Command Window only.', summary_log_file);
    dest_list = 1;
else
    dest_list = [1, fid_sum];
end

n_sub_half = 10; % 10 clean trials for Half A (basis) and 10 for Half B (projection)

for band_idx = 1:length(bands_to_process)
    band_name = bands_to_process{band_idx};
    current_data = data_all_conds.(band_name); 
    num_subjects = length(current_data.(cfg.ref_cond));
    
    for dest = dest_list
        fprintf(dest, '\n\n==================================================================================\n');
        fprintf(dest, '   BAND: %s | UNIVERSAL BASIS (m = %d VALIDATED DIRECTIONS)\n', upper(band_name), cfg.m);
        fprintf(dest, '==================================================================================\n');
    end
    
    group_dir = fullfile(output_path, 'group', band_name, 'Universal_5D');
    if ~exist(group_dir, 'dir'), mkdir(group_dir); end

    Xi_all_subjs = nan(cfg.num_ch, cfg.m, num_subjects);
    lam_frac_all = nan(cfg.m, num_subjects);
    gap_all      = nan(cfg.m - 1, num_subjects);
    lam_full_all = nan(cfg.num_ch, num_subjects); 
    
    for s = valid_retained_subjs
        rng(cfg.seed + s);
        
        % 1. Concatenate baseline-corrected Half-A ERPs (10 trials each) across all 5 conditions
        X_concat_A = [];
        for c_p = 1:length(cfg.all_conds)
            c_pname = cfg.all_conds{c_p};
            if ~isfield(current_data, c_pname), continue; end
            raw_c   = current_data.(c_pname){s};
            valid_c = find(~isnan(squeeze(raw_c(1, 1, :))));
            perm_c  = valid_c(randperm(length(valid_c)));
            n_h     = min(n_sub_half, floor(length(valid_c) / 2));
            
            erp_cA  = mean(raw_c(:, :, perm_c(1:n_h)), 3, 'omitnan');
            b_cA    = mean(erp_cA(:, iBase), 2, 'omitnan');
            X_cA_bc = erp_cA(:, iWin) - b_cA;
            X_concat_A = [X_concat_A, X_cA_bc]; %#ok<AGROW>
        end
        
        % 2. Centered dynamic covariance across the concatenated 5-condition matrix
        C_A_cent = covCentered(X_concat_A, zeros(cfg.num_ch, 1));
        
        [V, D] = eig(C_A_cent, 'vector');
        [lam, ord] = sort(D, 'descend');
        V = V(:, ord);
        
        % Retain exactly m = 5 orthogonal spatial modes
        Xi_all_subjs(:, :, s) = V(:, 1:cfg.m);
        lam_frac_all(:, s)    = lam(1:cfg.m) / sum(lam);
        gap_all(:, s)         = (lam(1:cfg.m-1) - lam(2:cfg.m)) ./ lam(1:cfg.m-1);
        lam_full_all(:, s)    = lam / sum(lam); 
    end
    
    aggregate_spectra.(band_name) = mean(cumsum(lam_full_all(:, valid_retained_subjs), 1), 2, 'omitnan') * 100;
    
    % Sign-align eigenvectors across participants against Subject 1
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
    
    % ------------------------------------------------------------------
    %% DELIVERABLE B: SPATIAL TOPOGRAPHIES (1x5 SUBPLOTS PER SUBJECT)
    % ------------------------------------------------------------------
    subjs_topo = valid_retained_subjs(1:min(4, length(valid_retained_subjs)));
    n_rows_B   = length(subjs_topo);
    
    figB = figure('Position', [100, 100, 220 * cfg.m + 120, 210 * n_rows_B], 'Color', 'w', 'Visible', 'off');
    for r_idx = 1:n_rows_B
        s = subjs_topo(r_idx);
        for k = 1:cfg.m
            ax_b = subplot(n_rows_B, cfg.m, (r_idx - 1) * cfg.m + k);
            pos_b = get(ax_b, 'Position');
            pos_b(1) = pos_b(1) * 0.92;
            pos_b(3) = pos_b(3) * 0.92;
            set(ax_b, 'Position', pos_b);
            
            topoplot(Xi_aligned(:, k, s), EEG.chanlocs, 'maplimits', [-0.4, 0.4], ...
                'electrodes', 'on', 'style', 'map', 'shading', 'interp');
            colormap(ax_b, 'jet'); clim(ax_b, [-0.4, 0.4]);
            
            if r_idx == 1
                title(sprintf('\\xi_%d', k), 'FontSize', 26, 'FontWeight', 'bold');
            end
            if k == 1
                text(-0.75, 0, sprintf('Subj %d', s), 'Rotation', 90, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 22, 'FontWeight', 'bold');
            end
        end
    end
    
    cb_B = colorbar('Position', [0.91, 0.15, 0.018, 0.72]);
    cb_B.Limits = [-0.4, 0.4]; cb_B.Ticks = -0.4 : 0.2 : 0.4;
    cb_B.Label.String = 'Spatial Weight (a.u.)'; cb_B.Label.FontSize = 22;
    cb_B.Label.FontWeight = 'bold'; cb_B.FontSize = 18;
    
    save_file_B = fullfile(group_dir, 'Deliverable_B_Topos');
    try saveas(figB, [save_file_B, '.png']); saveas(figB, [save_file_B, '.svg']); catch, end
    close(figB);

    % ------------------------------------------------------------------
    %% DELIVERABLE C: PAIRWISE COSINE SIMILARITY HEATMAPS (1x5 GRID)
    % ------------------------------------------------------------------
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
    end
    
    figC_sim = figure('Position', [100, 100, 310 * cfg.m + 120, 420], 'Color', 'w', 'Visible', 'off');
    tiledlayout(1, cfg.m, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for k = 1:cfg.m
        nexttile;
        imagesc(subj_subj_sim(:, :, k));
        colormap('jet'); clim([0 1]); axis square;
        set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'LineWidth', 1.2);
        title(sprintf('\\xi_%d', k), 'FontSize', 28, 'FontWeight', 'bold');
        
        if k == 1
            ylabel('Subject Index', 'FontSize', 24, 'FontWeight', 'bold');
        else
            yticklabels({});
        end
        xlabel('Subject Index', 'FontSize', 24, 'FontWeight', 'bold');
        xticks([1, n_ret]); yticks([1, n_ret]);
    end
    
    cb_sim = colorbar; cb_sim.Layout.Tile = 'east'; cb_sim.FontSize = 18;
    cb_sim.Label.String = 'Absolute Cosine Similarity'; cb_sim.Label.FontSize = 22; cb_sim.Label.FontWeight = 'bold';
    sgtitle(sprintf('Cross-Subject Component Consistency (%s Universal Basis, n=%d)', band_name, n_ret), ...
        'FontSize', 26, 'FontWeight', 'bold');
    
    save(fullfile(group_dir, 'subj_subj_similarity.mat'), 'subj_subj_sim', 'valid_retained_subjs', 'Xi_aligned');
    saveas(figC_sim, fullfile(group_dir, 'Deliverable_C_SubjSubj_Sim.png')); 
    close(figC_sim);

    % ------------------------------------------------------------------
    %% DELIVERABLE C (AUDIT): K-MEANS CLUSTERING (k=2 vs. k=3 ON DIRECTIONS 1 & 2)
    % ------------------------------------------------------------------
    X_feat = zeros(n_ret, 2 * cfg.num_ch);
    for i = 1:n_ret
        s = valid_retained_subjs(i);
        X_feat(i, :) = [Xi_aligned(:, 1, s)', Xi_aligned(:, 2, s)'];
    end
    
    sim_xi12 = 0.5 * (subj_subj_sim(:, :, 1) + subj_subj_sim(:, :, 2));
    
    rng(42);
    [idx_k2, ~] = kmeans(X_feat, 2, 'Replicates', 50, 'Distance', 'sqeuclidean');
    if sum(idx_k2 == 1) < sum(idx_k2 == 2), idx_k2 = 3 - idx_k2; end
    sil_k2 = mean(silhouette(X_feat, idx_k2));
    
    rng(42);
    [idx_k3, ~] = kmeans(X_feat, 3, 'Replicates', 50, 'Distance', 'sqeuclidean');
    [~, count_ord] = sort(histcounts(idx_k3, 1:4), 'descend');
    remap = zeros(3, 1); for r_i = 1:3, remap(count_ord(r_i)) = r_i; end
    idx_k3 = remap(idx_k3);
    sil_k3 = mean(silhouette(X_feat, idx_k3));
    
    % Permutation test for k=2 modular separation
    same_m2 = (idx_k2 == idx_k2') & ~eye(n_ret);
    diff_m2 = (idx_k2 ~= idx_k2');
    obs_diff_k2 = mean(sim_xi12(same_m2)) - mean(sim_xi12(diff_m2));
    null_diff_k2 = zeros(1000, 1);
    for p_c = 1:1000
        shuf_c2 = idx_k2(randperm(n_ret));
        s_m = (shuf_c2 == shuf_c2') & ~eye(n_ret);
        d_m = (shuf_c2 ~= shuf_c2');
        null_diff_k2(p_c) = mean(sim_xi12(s_m)) - mean(sim_xi12(d_m));
    end
    p_k2_perm = (1 + sum(null_diff_k2 >= obs_diff_k2)) / 1001;
    
    for dest = dest_list
        fprintf(dest, '\n--- K-MEANS SUBJECT CLUSTERING (DIRECTIONS 1 & 2) [%s] ---\n', upper(band_name));
        fprintf(dest, '  * K=2 Solution (Silhouette = %.3f, Permutation p = %.4f):\n', sil_k2, p_k2_perm);
        fprintf(dest, '      Cluster 1 (n=%d): Subjects %s\n', sum(idx_k2==1), mat2str(valid_retained_subjs(idx_k2==1)));
        fprintf(dest, '      Cluster 2 (n=%d): Subjects %s\n', sum(idx_k2==2), mat2str(valid_retained_subjs(idx_k2==2)));
        fprintf(dest, '  * K=3 Solution (Silhouette = %.3f):\n', sil_k3);
        fprintf(dest, '      Cluster 1 (n=%d): Subjects %s\n', sum(idx_k3==1), mat2str(valid_retained_subjs(idx_k3==1)));
        fprintf(dest, '      Cluster 2 (n=%d): Subjects %s\n', sum(idx_k3==2), mat2str(valid_retained_subjs(idx_k3==2)));
        fprintf(dest, '      Cluster 3 (n=%d): Subjects %s\n', sum(idx_k3==3), mat2str(valid_retained_subjs(idx_k3==3)));
    end
    
    figKMeans = figure('Position', [100, 100, 1550, 480], 'Color', 'w', 'Visible', 'off');
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    [~, ord_k2] = sort(idx_k2, 'ascend'); subjs_k2 = valid_retained_subjs(ord_k2);
    nexttile; hold on; imagesc(sim_xi12(ord_k2, ord_k2)); colormap('jet'); clim([0 1]); axis square;
    set(gca, 'YDir', 'reverse', 'FontSize', 12, 'FontWeight', 'bold');
    xline(sum(idx_k2==1) + 0.5, 'w--', 'LineWidth', 2.5); yline(sum(idx_k2==1) + 0.5, 'w--', 'LineWidth', 2.5);
    xticks(1:n_ret); xticklabels(arrayfun(@(x) sprintf('S%d', x), subjs_k2, 'UniformOutput', false)); xtickangle(45);
    yticks(1:n_ret); yticklabels(arrayfun(@(x) sprintf('S%d', x), subjs_k2, 'UniformOutput', false));
    title(sprintf('K=2 Clustered Sim (Sil = %.2f)', sil_k2), 'FontSize', 16, 'FontWeight', 'bold');
    
    [~, ord_k3] = sort(idx_k3, 'ascend'); subjs_k3 = valid_retained_subjs(ord_k3);
    nexttile; hold on; imagesc(sim_xi12(ord_k3, ord_k3)); colormap('jet'); clim([0 1]); axis square;
    set(gca, 'YDir', 'reverse', 'FontSize', 12, 'FontWeight', 'bold');
    xline(sum(idx_k3==1) + 0.5, 'w--', 'LineWidth', 2.5); yline(sum(idx_k3==1) + 0.5, 'w--', 'LineWidth', 2.5);
    xline(sum(idx_k3<=2) + 0.5, 'w--', 'LineWidth', 2.5); yline(sum(idx_k3<=2) + 0.5, 'w--', 'LineWidth', 2.5);
    xticks(1:n_ret); xticklabels(arrayfun(@(x) sprintf('S%d', x), subjs_k3, 'UniformOutput', false)); xtickangle(45);
    yticks(1:n_ret); yticklabels(arrayfun(@(x) sprintf('S%d', x), subjs_k3, 'UniformOutput', false));
    title(sprintf('K=3 Clustered Sim (Sil = %.2f)', sil_k3), 'FontSize', 16, 'FontWeight', 'bold');
    
    [~, score_2d] = pca(X_feat, 'NumComponents', 2);
    nexttile; hold on; box on; grid on; set(gca, 'FontSize', 14, 'FontWeight', 'bold');
    c_k2 = {[0.12 0.47 0.71], [0.84 0.15 0.16]};
    for i = 1:n_ret
        s_id = valid_retained_subjs(i);
        scatter(score_2d(i, 1), score_2d(i, 2), 140, 'filled', ...
            'MarkerFaceColor', c_k2{idx_k2(i)}, 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
        text(score_2d(i, 1) + 0.02, score_2d(i, 2) + 0.02, sprintf('S%d', s_id), ...
            'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.15 0.15 0.15]);
    end
    xlabel('Spatial Component 1', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Spatial Component 2', 'FontSize', 14, 'FontWeight', 'bold');
    title('2D Embedding (\xi_1 & \xi_2)', 'FontSize', 16, 'FontWeight', 'bold');
    
    cb_km = colorbar; cb_km.Layout.Tile = 'east'; cb_km.Label.String = 'Absolute Cosine Similarity'; cb_km.FontSize = 14;
    sgtitle(sprintf('Universal Basis Subject Clustering: K-Means (k=2 vs. k=3, %s)', band_name), ...
        'FontSize', 20, 'FontWeight', 'bold');
    saveas(figKMeans, fullfile(group_dir, 'Deliverable_C_KMeans_Clustering_K2_K3.png'));
    close(figKMeans);

    % ------------------------------------------------------------------
    %% 5-DIMENSIONAL SPLIT-HALF RESAMPLING & CONDITION PROJECTIONS
    % ------------------------------------------------------------------
    G_splits     = nan(length(cfg.all_conds), cfg.nRep, num_subjects);
    log_r_splits = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
    lat_splits   = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
    sust_splits  = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
    traj_splits  = nan(length(cfg.all_conds), cfg.m, num_windows, cfg.nRep, num_subjects);

    for subj = valid_retained_subjs
        subj_dir = fullfile(output_path, sprintf('subj%02d', subj), band_name, 'Universal_5D');
        if ~exist(subj_dir, 'dir'), mkdir(subj_dir); end
        save(fullfile(subj_dir, 'cfg.mat'), 'cfg');
        
        Xi_stable = Xi_aligned(:, :, subj);
        
        for rep = 1:cfg.nRep
            rng(cfg.seed + rep + subj*1000); 
            
            halfA_idx      = struct();
            halfB_idx      = struct();
            X_concat_rep_A = [];
            
            for c_p = 1:length(cfg.all_conds)
                cond_p = cfg.all_conds{c_p};
                if ~isfield(current_data, cond_p), continue; end
                raw_cp   = current_data.(cond_p){subj};
                valid_cp = find(~isnan(squeeze(raw_cp(1, 1, :))));
                perm_cp  = valid_cp(randperm(length(valid_cp)));
                
                n_sub = min(n_sub_half, floor(length(valid_cp) / 2));
                halfA_idx.(cond_p) = perm_cp(1:n_sub);
                halfB_idx.(cond_p) = perm_cp(n_sub + (1:n_sub));
                
                erp_cpA  = mean(raw_cp(:, :, halfA_idx.(cond_p)), 3, 'omitnan');
                b_cpA    = mean(erp_cpA(:, iBase), 2, 'omitnan');
                X_cpA_bc = erp_cpA(:, iWin) - b_cpA;
                X_concat_rep_A = [X_concat_rep_A, X_cpA_bc]; %#ok<AGROW>
            end
            
            C_A_cent = covCentered(X_concat_rep_A, zeros(cfg.num_ch, 1));
            [V, D]   = eig(C_A_cent, 'vector');
            [lam, ord] = sort(D, 'descend');
            V = V(:, ord);
            Xi_rep = V(:, 1:cfg.m); % 5 directions
            lam_frac_rep = lam(1:cfg.m) / sum(lam);
            gap_rep = (lam(1:cfg.m-1) - lam(2:cfg.m)) ./ lam(1:cfg.m-1);
            
            for k = 1:cfg.m
                if dot(Xi_rep(:, k), Xi_stable(:, k)) < 0
                    Xi_rep(:, k) = -Xi_rep(:, k);
                end
            end
            
            basis_vars = struct('Xi', Xi_rep, 'lam', lam, 'lam_frac', lam_frac_rep, 'gap', gap_rep);
            
            % Held-out Half B of Cued (P1) computes reference ceiling across 5 dimensions
            raw_ref = current_data.(cfg.ref_cond){subj};
            iA      = halfA_idx.(cfg.ref_cond);
            iB      = halfB_idx.(cfg.ref_cond);
            erpB    = mean(raw_ref(:, :, iB), 3, 'omitnan');
            bB      = mean(erpB(:, iBase), 2, 'omitnan');
            XB      = erpB(:, iWin) - bB;
            
            v_ceil = sum((Xi_rep' * XB).^2, 2) / T_win; % [5 x 1]
            V_ceil = sum(v_ceil);                       % Scalar 5D ceiling
            
            proj_B = Xi_rep' * XB; 
            [~, max_idx_B] = max(abs(proj_B), [], 2);
            t_lat_ceil = t_eval(max_idx_B); 
            
            split_vars = struct('iA', iA, 'iB', iB);
            proj_vars  = struct('v_ceil', v_ceil);
            
            for c = 1:length(cfg.all_conds)
                cond = cfg.all_conds{c};
                if ~isfield(current_data, cond), continue; end
                
                raw_cond = current_data.(cond){subj};
                idx_k    = halfB_idx.(cond); % Clean held-out Half B
                erp_k    = mean(raw_cond(:, :, idx_k), 3, 'omitnan');
                bk       = mean(erp_k(:, iBase), 2, 'omitnan');
                Xk       = erp_k(:, iWin) - bk;
                
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

    save(fullfile(group_dir, 'results.mat'), 'G_splits', 'log_r_splits', 'lat_splits', 'sust_splits', 'traj_splits', 'valid_retained_subjs');

   % =====================================================================
   %% CHECKPOINT 1 (B): PERMUTATION TEST FIGURE WITH NON-OVERLAPPING LABELS
    % =====================================================================
    figC_perm = figure('Position', [150, 150, 1550, 700], 'Color', 'w', ...
        'Name', sprintf('Permutation Test: %s', band_name), 'Visible', 'off');
    t_perm = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'normal');
    title(t_perm, sprintf('Cross-Subject Consistency: Null vs True (%s Universal Basis, n=%d)', ...
        band_name, n_ret), 'FontSize', 24, 'FontWeight', 'bold');
    
    for k = 1:cfg.m
        nexttile; hold on;
        set(gca, 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 1.2, 'TickDir', 'out');
        
        % 1. Null distribution histogram across 1,000 electrode permutations
        histogram(null_mean_sim(k, :), 30, 'Normalization', 'probability', ...
            'FaceColor', [0.72 0.72 0.72], 'EdgeColor', 'none');
        
        % 2. 95% One-tailed Null Threshold: Positioned strictly at TOP, LEFT of dashed line
        xline(thresh_95(k), 'k--', ...
            'LineWidth', 1.8, ...
            'Label', '95% Thresh', ...
            'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'top', ...
            'LabelHorizontalAlignment', 'left', ...
            'FontSize', 15);
        
        % 3. True Observed Mean Similarity: Shifted DOWN to MIDDLE, LEFT of solid line
        if p_vals(k) < 0.05
            line_color = [0.82 0.12 0.12]; % Crimson Red for statistically significant modes
        else
            line_color = [0.25 0.25 0.25]; % Charcoal Gray for non-significant modes
        end
        
        xline(true_mean_sim(k), '-', ...
            'Color', line_color, ...
            'LineWidth', 2.6, ...
            'Label', 'True Sim', ...
            'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'middle', ...
            'LabelHorizontalAlignment', 'left', ...
            'FontSize', 15);
        
        title(sprintf('\\xi_%d (p = %.3f)', k, p_vals(k)), 'FontSize', 20, 'FontWeight', 'bold');
        
        % Axis typography
        if k >= 4
            xlabel('Mean Abs. Cosine Similarity', 'FontSize', 17, 'FontWeight', 'bold');
        end
        if ismember(k, [1, 4])
            ylabel('Probability', 'FontSize', 18, 'FontWeight', 'bold');
        end
        
        % Headroom ceiling preventing right-edge clipping
        xlim([0, max([max(null_mean_sim(k, :)), true_mean_sim(k)]) * 1.18]);
        grid on; set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
        box on;
    end
    
    save_file_perm = fullfile(group_dir, 'Deliverable_C_Permutation_Test.png');
    try
        exportgraphics(figC_perm, save_file_perm, 'Resolution', 600);
    catch
        saveas(figC_perm, save_file_perm);
    end
    close(figC_perm);

    % ------------------------------------------------------------------
    %% CHECKPOINT 2, 3a & 3b: ALL RETAINED SUBJECTS TABLES & FIGURES
    % ------------------------------------------------------------------
    subj_mean_log_r = squeeze(mean(log_r_splits, 3, 'omitnan')); 
    subj_mean_traj  = squeeze(mean(traj_splits, 4, 'omitnan'));  
    
    all_subjs    = valid_retained_subjs;                     
    subj_groups  = {all_subjs};
    cohort_Vital = {sprintf('Average over %d retained subjects (ex. S3, S14)', length(all_subjs))};
    num_groups   = 1;
    
    time_mask_zoom = (t_centers >= -0.1) & (t_centers <= 1.0);
    time_mask_full = true(size(t_centers)); 
    
    global_clims = zeros(cfg.m, 1);
    all_log_r    = [];
    for k = 1:cfg.m
        max_c = 0;
        for g_idx = 1:num_groups
            grp_s = subj_groups{g_idx};
            if isempty(grp_s), continue; end
            Xi_grp_c = mean(Xi_aligned(:, :, grp_s), 3, 'omitnan');
            Xi_grp_c = Xi_grp_c ./ (vecnorm(Xi_grp_c) + eps);
            for c_idx = 1:length(cfg.all_conds)
                cond = cfg.all_conds{c_idx};
                erp_grp = mean(cat(3, current_data.(cond){grp_s}), 3, 'omitnan');
                max_c = max(max_c, max(abs(Xi_grp_c(:, k)' * erp_grp)) * 0.5);
                if k == 1
                    g_lr_m = mean(subj_mean_log_r(c_idx, :, grp_s), 3, 'omitnan');
                    g_lr_s = std(subj_mean_log_r(c_idx, :, grp_s), 0, 3, 'omitnan') ./ sqrt(length(grp_s));
                    all_log_r = [all_log_r; g_lr_m(:) + g_lr_s(:); g_lr_m(:) - g_lr_s(:)];
                end
            end
        end
        if max_c <= 0 || isnan(max_c), max_c = 10; end
        global_clims(k) = max_c;
    end
    if isempty(all_log_r), all_log_r = [-1.5, 1.5]; end
    cp2_y_min = min(-1.5, floor(min(all_log_r) * 1.15));
    cp2_y_max = max( 1.5, ceil(max(all_log_r) * 1.15));

    for g = 1:num_groups
        curr_subjs = subj_groups{g};
        curr_n = length(curr_subjs);
        g_name = subj_group_names{g};
        cohort_disp = strrep(g_name, '_', ' ');
        if curr_n == 0, continue; end
        
        grp_mean_log_r = mean(subj_mean_log_r(:, :, curr_subjs), 3, 'omitnan');      
        grp_se_log_r   = std(subj_mean_log_r(:, :, curr_subjs), 0, 3, 'omitnan') ./ sqrt(curr_n);
        
        aggregate_log_r.(band_name).cohort(g).mean = grp_mean_log_r;
        aggregate_log_r.(band_name).cohort(g).se   = grp_se_log_r;
        aggregate_log_r.(band_name).cohort(g).subj_data = subj_mean_log_r(:, :, curr_subjs);

        subj_mean_G = squeeze(mean(G_splits(:, :, curr_subjs), 2, 'omitnan')); 
        if curr_n == 1, subj_mean_G = subj_mean_G(:); end
        grp_mean_G = mean(subj_mean_G, 2, 'omitnan');
        aggregate_G.(band_name).cohort(g).mean = grp_mean_G;
        aggregate_G.(band_name).cohort(g).se   = std(subj_mean_G, 0, 2, 'omitnan') ./ sqrt(curr_n);
        aggregate_G.(band_name).cohort(g).subj_data = subj_mean_G;
        
        subj_mean_lat = squeeze(mean(mean(lat_splits(:, :, :, curr_subjs), 3, 'omitnan'), 2, 'omitnan')); 
        if curr_n == 1, subj_mean_lat = subj_mean_lat(:); end
        lat_dir_mean = squeeze(mean(lat_splits(:, :, :, curr_subjs), 3, 'omitnan'));
        if curr_n == 1, lat_dir_mean = reshape(lat_dir_mean, [length(cfg.all_conds), cfg.m, 1]); end
        grp_mean_lat_dir = mean(lat_dir_mean, 3, 'omitnan');
        
        aggregate_lat.(band_name).cohort(g).mean = grp_mean_lat_dir;
        aggregate_lat.(band_name).cohort(g).se   = std(lat_dir_mean, 0, 3, 'omitnan') ./ sqrt(curr_n);
        aggregate_lat.(band_name).cohort(g).subj_data = subj_mean_lat; 
        
        sust_dir_mean = squeeze(mean(sust_splits(:, :, :, curr_subjs), 3, 'omitnan')) * 100;
        if curr_n == 1, sust_dir_mean = reshape(sust_dir_mean, [length(cfg.all_conds), cfg.m, 1]); end
        grp_mean_sust_dir = mean(sust_dir_mean, 3, 'omitnan');
        
        subj_mean_sust = squeeze(mean(mean(sust_splits(:, :, :, curr_subjs), 3, 'omitnan'), 2, 'omitnan')) * 100;
        if curr_n == 1, subj_mean_sust = subj_mean_sust(:); end
        aggregate_sust.(band_name).cohort(g).mean = mean(subj_mean_sust, 2, 'omitnan');
        aggregate_sust.(band_name).cohort(g).se   = std(subj_mean_sust, 0, 2, 'omitnan') ./ sqrt(curr_n);
        aggregate_sust.(band_name).cohort(g).subj_data = subj_mean_sust; 
        
        % --- COMMAND WINDOW & LOG FILE SUMMARY TABLE PRINTOUT ---
        sig_components = find(p_vals < 0.05)';
        if isempty(sig_components), sig_components = 1:cfg.m; end
        
        for dest = dest_list
            fprintf(dest, '\n======================================================================\n');
            fprintf(dest, 'SUMMARY TABLE: %s REFERENCE (BAND: %s) [%s]\n', ...
                cfg.ref_cond, upper(band_name), cohort_Vital{g});
            fprintf(dest, '======================================================================\n');
            for k_idx = 1:length(sig_components)
                k = sig_components(k_idx);
                fprintf(dest, '\n--- SIGNIFICANT DIRECTION: \\xi_%d (Global p = %.3f) ---\n', k, p_vals(k));
                fprintf(dest, '----------------------------------------------------------------------\n');
                fprintf(dest, '%-10s | %-8s | %-8s | %-22s | %-18s\n', ...
                    'Condition', 'Gain (G)', sprintf('log(r_%d)', k), ...
                    sprintf('\\Delta Latency (\\xi_%d)', k), sprintf('Sustained %% (\\xi_%d)', k));
                fprintf(dest, '----------------------------------------------------------------------\n');
                for c_idx = 1:length(cfg.all_conds)
                    c_str = cfg.all_conds{c_idx};
                    fprintf(dest, '%-10s | %6.3f   | %+6.3f   | %+6.1f ms              | %6.2f%%\n', ...
                        c_str, grp_mean_G(c_idx), grp_mean_log_r(c_idx, k), ...
                        grp_mean_lat_dir(c_idx, k), grp_mean_sust_dir(c_idx, k));
                end
            end
        end
        
        % -----------------------------------------------------------------
        % CHECKPOINT 2: SUBSPACE REDISTRIBUTION PER BAND
        % -----------------------------------------------------------------
        figChk2 = figure('Position', [100, 100, 1400, 600], 'Visible', 'off');
        tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
        
        nexttile; hold on; set(gca, 'FontSize', 18);
        yline(0, 'k--', 'LineWidth', 1.8, 'DisplayName', sprintf('Ceiling (%s B)', cfg.ref_cond));
        colors_A = {[0.850 0.325 0.098], [0.929 0.694 0.125], [0.494 0.184 0.556]};
        for i = 1:length(cfg.group_A)
            c_idx = find(strcmp(cfg.all_conds, cfg.group_A{i}));
            if isempty(c_idx), continue; end
            col_idx = mod(i-1, length(colors_A)) + 1;
            errorbar(1:cfg.m, grp_mean_log_r(c_idx, :), grp_se_log_r(c_idx, :), '-s', ...
                'Color', colors_A{col_idx}, 'MarkerFaceColor', colors_A{col_idx}, ...
                'LineWidth', 2.5, 'MarkerSize', 8, 'DisplayName', clean_name(cfg.group_A{i}));
        end
        grid on; xlim([0.8, cfg.m + 0.2]); xticks(1:cfg.m); ylim([cp2_y_min, cp2_y_max]);
        xlabel('Spatial Direction Index', 'FontSize', 20); ylabel('log(r_i) \pm SEM', 'FontSize', 20);
        title(sprintf('Stimulus Delivered (%s)', cohort_disp), 'FontSize', 20, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 14);
        
        nexttile; hold on; set(gca, 'FontSize', 18);
        yline(0, 'k--', 'LineWidth', 1.8, 'DisplayName', sprintf('Ceiling (%s B)', cfg.ref_cond));
        colors_B = {[0.466 0.674 0.188], [0.301 0.745 0.933], [0.494 0.184 0.556]};
        grpB_sub = intersect({'P2_2000', 'P3_missing'}, cfg.group_B, 'stable');
        for i = 1:length(grpB_sub)
            c_idx = find(strcmp(cfg.all_conds, grpB_sub{i}));
            if isempty(c_idx), continue; end
            col_idx = mod(i-1, length(colors_B)) + 1;
            errorbar(1:cfg.m, grp_mean_log_r(c_idx, :), grp_se_log_r(c_idx, :), '-s', ...
                'Color', colors_B{col_idx}, 'MarkerFaceColor', colors_B{col_idx}, ...
                'LineWidth', 2.5, 'MarkerSize', 8, 'DisplayName', clean_name(grpB_sub{i}));
        end
        grid on; xlim([0.8, cfg.m + 0.2]); xticks(1:cfg.m); ylim([cp2_y_min, cp2_y_max]);
        xlabel('Spatial Direction Index', 'FontSize', 20); ylabel('log(r_i) \pm SEM', 'FontSize', 20);
        title(sprintf('Stimulus Omitted / Delayed (%s)', cohort_disp), 'FontSize', 20, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 14);
        sgtitle(sprintf('Subspace Redistribution: %s %s', band_name, cfg.ref_cond), 'FontSize', 24, 'FontWeight', 'bold');
        
        save_base_chk2 = fullfile(group_dir, sprintf('Checkpoint2_Redistribution_%s', g_name));
        try saveas(figChk2, [save_base_chk2, '.svg']); saveas(figChk2, [save_base_chk2, '.png']); catch, end
        close(figChk2);

        % -----------------------------------------------------------------
        % CHECKPOINT 3a: TRAJECTORIES (FULL-WINDOW Y-LIMITS + 15% HEADROOM)
        % -----------------------------------------------------------------
        grp_mean_traj = squeeze(mean(subj_mean_traj(:, :, :, curr_subjs), 4, 'omitnan')); 
        grp_sem_traj  = squeeze(std(subj_mean_traj(:, :, :, curr_subjs), 0, 4, 'omitnan')) ./ sqrt(curr_n);
        
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
            cp3a_max_y = max(row_upper(:)) * 1.15; 
            if isempty(cp3a_max_y) || isnan(cp3a_max_y) || cp3a_max_y <= cp3a_min_y
                cp3a_min_y = 0; cp3a_max_y = 10;
            end
            
            % --- Left Column: Stimulus Delivered ---
            nexttile; hold on; set(gca, 'FontSize', 22);
            patch([cfg.win(1) cfg.win(2) cfg.win(2) cfg.win(1)], [cp3a_min_y cp3a_min_y cp3a_max_y cp3a_max_y], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
            xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
            
            if ~isempty(c_idx_ref)
                y_val_ref = reshape(grp_mean_traj(c_idx_ref, k, :), 1, []); y_sem_ref = reshape(grp_sem_traj(c_idx_ref, k, :), 1, []);
                fill([t_centers, fliplr(t_centers)], max(0, [y_val_ref + y_sem_ref, fliplr(y_val_ref - y_sem_ref)]), [0.5 0.5 0.5], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(t_centers, y_val_ref, 'k-', 'LineWidth', 2, 'DisplayName', sprintf('Reference (%s)', clean_name(cfg.ref_cond)));
            end
            
            for i = 1:length(cfg.group_A)
                c_idx = find(strcmp(cfg.all_conds, cfg.group_A{i}));
                if isempty(c_idx), continue; end
                y_val = reshape(grp_mean_traj(c_idx, k, :), 1, []); y_sem = reshape(grp_sem_traj(c_idx, k, :), 1, []);
                col_idx = mod(i-1, 3) + 1;
                fill([t_centers, fliplr(t_centers)], max(0, [y_val + y_sem, fliplr(y_val - y_sem)]), colors_A{col_idx}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(t_centers, y_val, '-', 'Color', colors_A{col_idx}, 'LineWidth', 2, 'DisplayName', clean_name(cfg.group_A{i}));
            end
            grid on; xlabel('Time (s)', 'FontSize', 20); ylabel('Projected Variance', 'FontSize', 20);
            title(sprintf('Stimulus Delivered: \\xi_%d Trajectory', k), 'FontSize', 24); 
            if k_idx == 1, lgdA = legend('Location', 'best'); lgdA.FontSize = 18; end
            xlim([t_centers(1), t_centers(end)]); ylim([cp3a_min_y, cp3a_max_y]);
            
            % --- Right Column: Stimulus Omitted / Delayed ---
            nexttile; hold on; set(gca, 'FontSize', 22);
            patch([cfg.win(1) cfg.win(2) cfg.win(2) cfg.win(1)], [cp3a_min_y cp3a_min_y cp3a_max_y cp3a_max_y], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
            xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
            
            if ~isempty(c_idx_ref)
                fill([t_centers, fliplr(t_centers)], max(0, [y_val_ref + y_sem_ref, fliplr(y_val_ref - y_sem_ref)]), [0.5 0.5 0.5], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(t_centers, y_val_ref, 'k-', 'LineWidth', 2, 'DisplayName', sprintf('Reference (%s)', clean_name(cfg.ref_cond)));
            end
            
            for i = 1:length(cfg.group_B)
                c_idx = find(strcmp(cfg.all_conds, cfg.group_B{i}));
                if isempty(c_idx), continue; end
                y_val = reshape(grp_mean_traj(c_idx, k, :), 1, []); y_sem = reshape(grp_sem_traj(c_idx, k, :), 1, []);
                col_idx = mod(i-1, length(colors_B)) + 1;
                fill([t_centers, fliplr(t_centers)], max(0, [y_val + y_sem, fliplr(y_val - y_sem)]), colors_B{col_idx}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(t_centers, y_val, '-', 'Color', colors_B{col_idx}, 'LineWidth', 2, 'DisplayName', clean_name(cfg.group_B{i}));
            end
            grid on; xlabel('Time (s)', 'FontSize', 20); ylabel('Projected Variance', 'FontSize', 20);
            title(sprintf('Stimulus Omitted / Delayed: \\xi_%d Trajectory', k), 'FontSize', 24); 
            if k_idx == 1, lgdB = legend('Location', 'best'); lgdB.FontSize = 18; end
            xlim([t_centers(1), t_centers(end)]); ylim([cp3a_min_y, cp3a_max_y]);
        end
        sgtitle(sprintf('Subspace Trajectories: %s (%s, n=%d)', band_name, cohort_disp, curr_n), 'FontSize', 24, 'FontWeight', 'bold');
        
        save_base_chk3 = fullfile(group_dir, sprintf('Checkpoint3a_Trajectories_%s', g_name));
        try saveas(figChk3, [save_base_chk3, '.svg']); saveas(figChk3, [save_base_chk3, '.png']); catch, end
        close(figChk3);

        % =========================================================================
        % CHECKPOINT 3b: SINGLE TOPOPLOT + TRAJECTORY FILMSTRIPS (ALL 6 DIRECTIONS)
        % =========================================================================
        orig_warn_state = warning('off', 'all'); 

        Xi_bar = mean(Xi_aligned(:, :, curr_subjs), 3, 'omitnan');
        Xi_bar = Xi_bar ./ vecnorm(Xi_bar);

        topo_win_s = [0.2, 0.8];
        topo_idx = find(time_s >= topo_win_s(1) & time_s <= topo_win_s(2));

        groupA_conds = {'P1', 'P2_500', 'P3_500'};
        groupB_conds = {'P1', 'P2_2000', 'P3_missing'};
        colors_A = {
            [0 0 0], ...
            [0.850 0.325 0.098], ...
            [0.929 0.694 0.125]
        };
        colors_B = {
            [0 0 0], ...
            [0.466 0.674 0.188], ...
            [0.301 0.745 0.933], ...
            [0.494 0.184 0.556]
        };

        for k = 1:cfg.m
            xi = Xi_bar(:, k);

            time_mask_zoom = (t_centers >= -0.1 & t_centers <= 1.0);
            zoom_vals_up = reshape( ...
                grp_mean_traj(:, k, time_mask_zoom) + ...
                grp_sem_traj(:, k, time_mask_zoom), [], 1);
            zoom_vals_dn = reshape( ...
                grp_mean_traj(:, k, time_mask_zoom) - ...
                grp_sem_traj(:, k, time_mask_zoom), [], 1);
            max_y = max(zoom_vals_up, [], 'omitnan') * 1.15;
            min_y = max(0, min(zoom_vals_dn, [], 'omitnan') * 0.90);
            if isempty(max_y) || ~isfinite(max_y) || ...
                    isempty(min_y) || ~isfinite(min_y) || max_y <= min_y
                min_y = 0;
                max_y = 10;
            end

            if isfield(current_data, cfg.ref_cond) && ...
                    ~isempty(current_data.(cfg.ref_cond))
                erp_cued = cat(3, current_data.(cfg.ref_cond){curr_subjs});
                avg_cued = mean(erp_cued, 3, 'omitnan');
                cued_topo_vector = xi * mean( ...
                    xi' * avg_cued(:, topo_idx), 2, 'omitnan');
                c_lim = [
                    min(cued_topo_vector(:), [], 'omitnan'), ...
                    max(cued_topo_vector(:), [], 'omitnan')
                ];
            else
                c_lim = [-1, 1];
            end

            if any(~isfinite(c_lim)) || diff(c_lim) <= eps(max(abs(c_lim)))
                c_lim = [-1, 1];
            else
                max_abs_c = max(abs(c_lim));
                c_lim = [-max_abs_c, max_abs_c];
            end

            % -----------------------------------------------------------------
            % FIGURE 1: STIMULUS DELIVERED
            % -----------------------------------------------------------------
            n_conds_A = length(groupA_conds);
            figA = figure( ...
                'Position', [50, 50, 900, 260 * n_conds_A], ...
                'Color', 'w', ...
                'Name', sprintf('Dir %d - Stimulus Delivered (%s)', k, cohort_disp), ...
                'Visible', 'off');
            tA = tiledlayout(n_conds_A, 3, ...
                'TileSpacing', 'compact', ...
                'Padding', 'normal');

            for c = 1:n_conds_A
                cond = groupA_conds{c};
                col = colors_A{c};
                c_name_clean = clean_name(cond);

                erp_all = cat(3, current_data.(cond){curr_subjs});
                avg_data = mean(erp_all, 3, 'omitnan');
                topo_vector = xi * mean( ...
                    xi' * avg_data(:, topo_idx), 2, 'omitnan');
                topo_vector = topo_vector(:);

                ax_topo = nexttile(tA, [1, 1]);
                topoplot(topo_vector, EEG.chanlocs, ...
                    'maplimits', c_lim, ...
                    'electrodes', 'on', ...
                    'style', 'map', ...
                    'shading', 'interp');
                
                xlim(ax_topo, [-0.58, 0.58]);
                ylim(ax_topo, [-0.58, 0.62]);
                colormap(ax_topo, 'jet');
                clim(ax_topo, c_lim);
                title(ax_topo, sprintf('%s (0.2-0.8s)', c_name_clean), ...
                    'FontSize', 12, 'FontWeight', 'bold');

                cb = colorbar(ax_topo, 'eastoutside');
                cb.Label.String = 'Amp (a.u.)';
                cb.Label.FontSize = 10;

                ax_traj = nexttile(tA, [1, 2]);
                hold(ax_traj, 'on');
                patch(ax_traj, ...
                    [cfg.win(1) cfg.win(2) cfg.win(2) cfg.win(1)], ...
                    [min_y min_y max_y max_y], ...
                    [0.9 0.9 0.9], ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.5, ...
                    'HandleVisibility', 'off');
                xline(ax_traj, 0, 'k--', ...
                    'LineWidth', 1, ...
                    'HandleVisibility', 'off');

                c_idx = find(strcmp(cfg.all_conds, cond));
                if ~isempty(c_idx)
                    y_val = reshape(grp_mean_traj(c_idx, k, :), 1, []);
                    y_sem = reshape(grp_sem_traj(c_idx, k, :), 1, []);
                    fill(ax_traj, ...
                        [t_centers, fliplr(t_centers)], ...
                        max(0, [y_val + y_sem, ...
                                fliplr(y_val - y_sem)]), ...
                        col, ...
                        'FaceAlpha', 0.25, ...
                        'EdgeColor', 'none', ...
                        'HandleVisibility', 'off');
                    plot(ax_traj, t_centers, y_val, '-', ...
                        'Color', col, ...
                        'LineWidth', 2.5);
                end

                xlim(ax_traj, [-0.1, 1.0]);
                ylim(ax_traj, [min_y, max_y]);
                ylabel(ax_traj, 'Proj. Var', ...
                    'FontSize', 11, 'FontWeight', 'bold');
                grid(ax_traj, 'on');
                title(ax_traj, sprintf('%s - Trajectory', c_name_clean), ...
                    'FontSize', 12, 'FontWeight', 'bold');
                if c == n_conds_A
                    xlabel(ax_traj, 'Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
                end
            end

            title(tA, ...
                sprintf('Stimulus Delivered: Direction (\\xi_%d) [%s | %s]', ...
                k, band_name, cohort_disp), ...
                'FontSize', 16, ...
                'FontWeight', 'bold');

            save_file_A = fullfile(group_dir, ...
                sprintf('Checkpoint3b_GroupA_Summary_Xi%d_%s.png', k, g_name));
            try
                exportgraphics(figA, save_file_A, 'Resolution', 300);
            catch
                try saveas(figA, save_file_A); catch, warning('Could not save %s.', save_file_A); end
            end
            close(figA);

            % -----------------------------------------------------------------
            % FIGURE 2: STIMULUS OMITTED / DELAYED
            % -----------------------------------------------------------------
            n_conds_B = length(groupB_conds);
            figB = figure( ...
                'Position', [100, 100, 900, 260 * n_conds_B], ...
                'Color', 'w', ...
                'Name', sprintf('Dir %d - Stimulus Omitted or Delayed (%s)', k, cohort_disp), ...
                'Visible', 'off');
            tB = tiledlayout(n_conds_B, 3, ...
                'TileSpacing', 'compact', ...
                'Padding', 'normal');

            for c = 1:n_conds_B
                cond = groupB_conds{c};
                col = colors_B{c};
                c_name_clean = clean_name(cond);

                erp_all = cat(3, current_data.(cond){curr_subjs});
                avg_data = mean(erp_all, 3, 'omitnan');
                topo_vector = xi * mean( ...
                    xi' * avg_data(:, topo_idx), 2, 'omitnan');
                topo_vector = topo_vector(:);

                ax_topo = nexttile(tB, [1, 1]);
                topoplot(topo_vector, EEG.chanlocs, ...
                    'maplimits', c_lim, ...
                    'electrodes', 'on', ...
                    'style', 'map', ...
                    'shading', 'interp');
                
                xlim(ax_topo, [-0.58, 0.58]);
                ylim(ax_topo, [-0.58, 0.62]);
                colormap(ax_topo, 'jet');
                clim(ax_topo, c_lim);
                title(ax_topo, sprintf('%s (0.2-0.8s)', c_name_clean), ...
                    'FontSize', 12, 'FontWeight', 'bold');

                cb = colorbar(ax_topo, 'eastoutside');
                cb.Label.String = 'Amp (a.u.)';
                cb.Label.FontSize = 10;

                ax_traj = nexttile(tB, [1, 2]);
                hold(ax_traj, 'on');
                patch(ax_traj, ...
                    [cfg.win(1) cfg.win(2) cfg.win(2) cfg.win(1)], ...
                    [min_y min_y max_y max_y], ...
                    [0.9 0.9 0.9], ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.5, ...
                    'HandleVisibility', 'off');
                xline(ax_traj, 0, 'k--', ...
                    'LineWidth', 1, ...
                    'HandleVisibility', 'off');

                c_idx = find(strcmp(cfg.all_conds, cond));
                if ~isempty(c_idx)
                    y_val = reshape(grp_mean_traj(c_idx, k, :), 1, []);
                    y_sem = reshape(grp_sem_traj(c_idx, k, :), 1, []);
                    fill(ax_traj, ...
                        [t_centers, fliplr(t_centers)], ...
                        max(0, [y_val + y_sem, ...
                                fliplr(y_val - y_sem)]), ...
                        col, ...
                        'FaceAlpha', 0.25, ...
                        'EdgeColor', 'none', ...
                        'HandleVisibility', 'off');
                    plot(ax_traj, t_centers, y_val, '-', ...
                        'Color', col, ...
                        'LineWidth', 2.5);
                end

                xlim(ax_traj, [-0.1, 1.0]);
                ylim(ax_traj, [min_y, max_y]);
                ylabel(ax_traj, 'Proj. Var', ...
                    'FontSize', 11, 'FontWeight', 'bold');
                grid(ax_traj, 'on');
                title(ax_traj, sprintf('%s - Trajectory', c_name_clean), ...
                    'FontSize', 12, 'FontWeight', 'bold');
                if c == n_conds_B
                    xlabel(ax_traj, 'Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
                end
            end

            title(tB, ...
                sprintf('Stimulus Omitted / Delayed: Direction (\\xi_%d) [%s | %s]', ...
                k, band_name, cohort_disp), ...
                'FontSize', 16, ...
                'FontWeight', 'bold');

            save_file_B = fullfile(group_dir, ...
                sprintf('Checkpoint3b_GroupB_Summary_Xi%d_%s.png', k, g_name));
            try
                exportgraphics(figB, save_file_B, 'Resolution', 300);
            catch
                try saveas(figB, save_file_B); catch, warning('Could not save %s.', save_file_B); end
            end
            close(figB);
        end

        warning(orig_warn_state);
    end
end

fclose(fid_sum);
fprintf('\nAll Summary Tables printed and saved to: %s\n', summary_log_file);

%% =========================================================================
%% FIGURE 1B: CHANNEL C3 ERP DYNAMICS + TOPOPLOTS (150, 400, 900 ms) + WILCOXON BARS
%% =========================================================================
fprintf('\nGenerating Figure 1B: Channel C3 ERPs with Topoplots (150, 400, 900 ms) & Wilcoxon Bars...\n');
orig_warn_1B = warning('off', 'all');

% -------------------------------------------------------------------------
% 1. USER CONFIGURATION & STATISTICAL PARAMETERS
% -------------------------------------------------------------------------
time_win       = [-0.200, 1.000];        % ERP time window [-0.2 s to 1.0 s]
topo_times     = [0.150, 0.400, 0.900];  % 3 topoplot latencies: 150, 400, 900 ms
topo_half_win  = 0.025;                  % Average +/- 25 ms around each latency for clean SNR
norm_au        = false;                  % false = Voltage (\muV); true = Normalized to P1 peak (a.u.)

% Statistical bar settings (Wilcoxon signed-rank against Cued P1)
alpha_stat     = 0.05;                   % Significance threshold (p < 0.05)
use_fdr        = false;                  % false = uncorrected + cluster duration; true = Benjamini-Hochberg FDR
min_sig_ms     = 20;                     % Minimum contiguous significant duration in ms
single_red_bar = false;                  % false = condition-colored bars; true = single crimson bar

if ~exist('chan_c3', 'var') || isempty(chan_c3)
    chan_c3 = 14;
    if exist('EEG', 'var') && isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs)
        match_idx = find(strcmpi({EEG.chanlocs.labels}, 'C3'), 1);
        if ~isempty(match_idx), chan_c3 = match_idx; end
    end
end

subjs_plot   = valid_retained_subjs;
n_subjs_plot = length(subjs_plot);
time_mask    = (time_s >= time_win(1)) & (time_s <= time_win(2));
t_eval       = time_s(time_mask);
n_t          = length(t_eval);

% -------------------------------------------------------------------------
% 2. COLOR PALETTE
% -------------------------------------------------------------------------
c_P1      = [0.07, 0.35, 0.37]; % Dark Teal
c_P2_500  = [0.47, 0.12, 0.08]; % Dark Maroon
c_P3_500  = [0.88, 0.57, 0.10]; % Amber / Gold
c_P2_2000 = [0.90, 0.59, 0.51]; % Coral / Salmon
c_P3_null = [0.94, 0.75, 0.47]; % Light Gold / Sand
c_sig_red = [0.78, 0.12, 0.12]; % Crimson Red

% -------------------------------------------------------------------------
% 3. EXTRACT C3 SUBJECT TRACES & 32-CHANNEL GRAND AVERAGES FOR TOPOPLOTS
% -------------------------------------------------------------------------
cond_list = {'P1', 'P2_500', 'P3_500', 'P2_2000', 'P3_missing'};
erp_stats = struct();

for c_i = 1:length(cond_list)
    c_name      = cond_list{c_i};
    subj_traces = zeros(n_t, n_subjs_plot);
    subj_32ch   = zeros(cfg.num_ch, length(time_s), n_subjs_plot);
    
    for s_i = 1:n_subjs_plot
        s       = subjs_plot(s_i);
        raw_mat = data_all_conds.Raw.(c_name){s}; % [32 x time x trials]
        
        erp_s = mean(raw_mat, 3, 'omitnan');
        subj_32ch(:, :, s_i) = erp_s;
        subj_traces(:, s_i)  = erp_s(chan_c3, time_mask)';
    end
    
    erp_stats.(c_name).subj_traces = subj_traces;
    erp_stats.(c_name).mean        = mean(subj_traces, 2, 'omitnan');
    erp_stats.(c_name).topo_32ch   = mean(subj_32ch, 3, 'omitnan');
end

if norm_au
    p1_peak = max(abs(erp_stats.P1.mean));
    if p1_peak <= 0, p1_peak = 1; end
    for c_i = 1:length(cond_list)
        c_name = cond_list{c_i};
        erp_stats.(c_name).mean        = erp_stats.(c_name).mean / p1_peak;
        erp_stats.(c_name).subj_traces = erp_stats.(c_name).subj_traces / p1_peak;
        erp_stats.(c_name).topo_32ch   = erp_stats.(c_name).topo_32ch / p1_peak;
    end
    y_label_str = 'EEG Amplitude (a.u.)';
    cb_unit_str = 'a.u.';
else
    y_label_str = 'EEG Amplitude (\muV)';
    cb_unit_str = '\muV';
end

% -------------------------------------------------------------------------
% 4. PAIRED WILCOXON SIGNED-RANK TESTS ACROSS TIME AT CHANNEL C3
% -------------------------------------------------------------------------
pairs_A  = {'P2_500', 'P3_500'};
pairs_B  = {'P2_2000', 'P3_missing'};
p_vals_A = ones(length(pairs_A), n_t);
p_vals_B = ones(length(pairs_B), n_t);

for p_i = 1:length(pairs_A)
    c_comp = pairs_A{p_i};
    for t_i = 1:n_t
        x_ref  = erp_stats.P1.subj_traces(t_i, :);
        x_comp = erp_stats.(c_comp).subj_traces(t_i, :);
        if any(~isnan(x_ref - x_comp))
            p_vals_A(p_i, t_i) = signrank(x_comp, x_ref);
        end
    end
end

for p_i = 1:length(pairs_B)
    c_comp = pairs_B{p_i};
    for t_i = 1:n_t
        x_ref  = erp_stats.P1.subj_traces(t_i, :);
        x_comp = erp_stats.(c_comp).subj_traces(t_i, :);
        if any(~isnan(x_ref - x_comp))
            p_vals_B(p_i, t_i) = signrank(x_comp, x_ref);
        end
    end
end

min_pts      = max(1, round((min_sig_ms / 1000) * cfg.fs));
get_sig_mask = @(p_vec) filter_contiguous_runs(apply_bh_fdr(p_vec, alpha_stat, use_fdr), min_pts);

sig_A = false(length(pairs_A), n_t);
sig_B = false(length(pairs_B), n_t);
for p_i = 1:length(pairs_A), sig_A(p_i, :) = get_sig_mask(p_vals_A(p_i, :)); end
for p_i = 1:length(pairs_B), sig_B(p_i, :) = get_sig_mask(p_vals_B(p_i, :)); end

% -------------------------------------------------------------------------
% 5. COMPUTE COORDINATED Y-LIMITS & TOPOPLOT COLOR LIMITS
% -------------------------------------------------------------------------
all_means = [];
for c_i = 1:length(cond_list)
    all_means = [all_means; erp_stats.(cond_list{c_i}).mean];
end
min_val = min(all_means);
max_val = max(all_means);
y_span  = max_val - min_val;
if y_span <= 0 || isnan(y_span), y_span = 2; end
ylim_shared = [min_val - 0.24 * y_span, max_val + 0.15 * y_span];

topo_map_A = (erp_stats.P1.topo_32ch + erp_stats.P2_500.topo_32ch + erp_stats.P3_500.topo_32ch) / 3;
topo_map_B = (erp_stats.P2_2000.topo_32ch + erp_stats.P3_missing.topo_32ch) / 2;

topo_vecs_A = zeros(cfg.num_ch, length(topo_times));
topo_vecs_B = zeros(cfg.num_ch, length(topo_times));
for k = 1:length(topo_times)
    idx_w = (time_s >= (topo_times(k) - topo_half_win)) & (time_s <= (topo_times(k) + topo_half_win));
    topo_vecs_A(:, k) = mean(topo_map_A(:, idx_w), 2, 'omitnan');
    topo_vecs_B(:, k) = mean(topo_map_B(:, idx_w), 2, 'omitnan');
end
max_topo = max(abs([topo_vecs_A(:); topo_vecs_B(:)])) * 0.90;
if isempty(max_topo) || max_topo <= 0 || isnan(max_topo), max_topo = 1; end
c_lim_topo = [-max_topo, max_topo];

% -------------------------------------------------------------------------
% 6. RENDER PUBLICATION FIGURE WITH BALANCED FONTS & ADJUSTED SPACING
% -------------------------------------------------------------------------
fig1B = figure('Position', [50, 50, 1500, 660], 'Color', 'w', ...
    'Name', 'Figure 1B: Channel C3 ERPs with Topoplots & Wilcoxon Bars');

% Adjusted subplot coordinates to widen the inter-panel gap
ax1_pos = [0.065, 0.13, 0.380, 0.54];
ax2_pos = [0.575, 0.13, 0.380, 0.54];

% ==================== SUBPLOT 1: STIMULUS DELIVERED ====================
ax1 = axes('Position', ax1_pos); hold(ax1, 'on');
% Decreased axis font size from 26 to 18 to prevent tick collisions and auto-rotation
set(ax1, 'FontSize', 18, 'FontName', 'Helvetica', 'LineWidth', 1.5, 'TickDir', 'out');
box(ax1, 'off');

yline(ax1, 0, 'Color', [0.82 0.82 0.82], 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(ax1, 0, 'k--', 'LineWidth', 1.2, 'Alpha', 0.6, 'HandleVisibility', 'off');

for k = 1:length(topo_times)
    tk = topo_times(k);
    plot(ax1, [tk, tk], [ylim_shared(1), ylim_shared(2) * 1.05], 'k:', ...
        'LineWidth', 1.2, 'Clipping', 'off', 'HandleVisibility', 'off');
end

plot(ax1, t_eval, erp_stats.P1.mean,     'Color', c_P1,     'LineWidth', 2.6, 'DisplayName', 'Cued (P1)');
plot(ax1, t_eval, erp_stats.P2_500.mean, 'Color', c_P2_500, 'LineWidth', 2.6, 'DisplayName', 'Unpred 500 (P2-500)');
plot(ax1, t_eval, erp_stats.P3_500.mean, 'Color', c_P3_500, 'LineWidth', 2.6, 'DisplayName', 'Rand 500 (P3-500)');

% Wilcoxon Significance Bars
y_bar_base = ylim_shared(1) + 0.04 * (ylim_shared(2) - ylim_shared(1));
y_bar_step = 0.045 * (ylim_shared(2) - ylim_shared(1));

if single_red_bar
    any_sig_A = any(sig_A, 1);
    draw_sig_bars(ax1, t_eval, any_sig_A, y_bar_base, c_sig_red, 6.0);
    if any(any_sig_A)
        text(ax1, 0.52, y_bar_base + y_bar_step*0.8, 'p < 0.05 (Wilcoxon)', ...
            'Color', c_sig_red, 'FontSize', 12, 'FontWeight', 'bold');
    end
else
    draw_sig_bars(ax1, t_eval, sig_A(1, :), y_bar_base + y_bar_step, c_P2_500, 5.0);
    draw_sig_bars(ax1, t_eval, sig_A(2, :), y_bar_base,              c_P3_500, 5.0);
    text(ax1, -0.18, y_bar_base + y_bar_step*0.5, 'p < 0.05 vs. P1', ...
        'Color', [0.25 0.25 0.25], 'FontSize', 11, 'FontWeight', 'bold');
end

text(ax1, 0.82, ylim_shared(2)*0.86, 'P1',     'Color', c_P1,     'FontSize', 19, 'FontWeight', 'bold');
text(ax1, 0.82, ylim_shared(2)*0.72, 'P2-500', 'Color', c_P2_500, 'FontSize', 19, 'FontWeight', 'bold');
text(ax1, 0.82, ylim_shared(2)*0.58, 'P3-500', 'Color', c_P3_500, 'FontSize', 19, 'FontWeight', 'bold');

xlim(ax1, time_win); ylim(ax1, ylim_shared);
xticks(ax1, -0.2 : 0.2 : 1.0);
xticklabels(ax1, arrayfun(@(x) sprintf('%d', round(x*1000)), -0.2:0.2:1.0, 'UniformOutput', false));
xtickangle(ax1, 0); % Enforce horizontal x-axis tick labels
xlabel(ax1, 'Time after expected stimulus (ms)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel(ax1, y_label_str, 'FontSize', 20, 'FontWeight', 'bold');

% ==================== SUBPLOT 2: STIMULUS OMITTED / DELAYED ====================
ax2 = axes('Position', ax2_pos); hold(ax2, 'on');
% Increased axis font size from 15 to 18 to match the left subplot perfectly
set(ax2, 'FontSize', 18, 'FontName', 'Helvetica', 'LineWidth', 1.5, 'TickDir', 'out');
box(ax2, 'off');

yline(ax2, 0, 'Color', [0.82 0.82 0.82], 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(ax2, 0, 'k--', 'LineWidth', 1.2, 'Alpha', 0.6, 'HandleVisibility', 'off');

for k = 1:length(topo_times)
    tk = topo_times(k);
    plot(ax2, [tk, tk], [ylim_shared(1), ylim_shared(2) * 1.05], 'k:', ...
        'LineWidth', 1.2, 'Clipping', 'off', 'HandleVisibility', 'off');
end

plot(ax2, t_eval, erp_stats.P1.mean,         'Color', c_P1,      'LineWidth', 2.6, 'DisplayName', 'Cued (P1)');
plot(ax2, t_eval, erp_stats.P2_2000.mean,    'Color', c_P2_2000, 'LineWidth', 2.6, 'DisplayName', 'Unpred 2000 (P2-2000)');
plot(ax2, t_eval, erp_stats.P3_missing.mean, 'Color', c_P3_null, 'LineWidth', 2.6, 'DisplayName', 'Rand Null (P3-null)');

if single_red_bar
    any_sig_B = any(sig_B, 1);
    draw_sig_bars(ax2, t_eval, any_sig_B, y_bar_base, c_sig_red, 6.0);
    if any(any_sig_B)
        text(ax2, 0.52, y_bar_base + y_bar_step*0.8, 'p < 0.05 (Wilcoxon)', ...
            'Color', c_sig_red, 'FontSize', 12, 'FontWeight', 'bold');
    end
else
    draw_sig_bars(ax2, t_eval, sig_B(1, :), y_bar_base + y_bar_step, c_P2_2000, 5.0);
    draw_sig_bars(ax2, t_eval, sig_B(2, :), y_bar_base,              c_P3_null, 5.0);
    text(ax2, -0.18, y_bar_base + y_bar_step*0.5, 'p < 0.05 vs. P1', ...
        'Color', [0.25 0.25 0.25], 'FontSize', 11, 'FontWeight', 'bold');
end

text(ax2, 0.80, ylim_shared(2)*0.86, 'P1',      'Color', c_P1,      'FontSize', 19, 'FontWeight', 'bold');
text(ax2, 0.80, ylim_shared(2)*0.72, 'P2-2000', 'Color', c_P2_2000, 'FontSize', 19, 'FontWeight', 'bold');
text(ax2, 0.80, ylim_shared(2)*0.58, 'P3-null', 'Color', c_P3_null, 'FontSize', 19, 'FontWeight', 'bold');

xlim(ax2, time_win); ylim(ax2, ylim_shared);
xticks(ax2, -0.2 : 0.2 : 1.0);
xticklabels(ax2, arrayfun(@(x) sprintf('%d', round(x*1000)), -0.2:0.2:1.0, 'UniformOutput', false));
xtickangle(ax2, 0); % Enforce horizontal x-axis tick labels
xlabel(ax2, 'Time after expected stimulus (ms)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel(ax2, y_label_str, 'FontSize', 20, 'FontWeight', 'bold');

% -------------------------------------------------------------------------
% 7. FLOATING TOPOPLOTS & NON-OVERLAPPING COLORBAR
% -------------------------------------------------------------------------
topo_w = 0.076;                           
topo_h = 0.175;                           
topo_y = ax1_pos(2) + ax1_pos(4) + 0.015; 
t_span = time_win(2) - time_win(1);

% Panel A Topoplots (150, 400, 900 ms)
for k = 1:length(topo_times)
    frac_x   = (topo_times(k) - time_win(1)) / t_span;
    center_x = ax1_pos(1) + frac_x * ax1_pos(3);
    ax_tA    = axes('Position', [center_x - topo_w/2, topo_y, topo_w, topo_h]); %#ok<LAXES>
    topoplot(topo_vecs_A(:, k), EEG.chanlocs, 'maplimits', c_lim_topo, ...
        'numcontour', 3, 'electrodes', 'on', 'style', 'map', 'shading', 'interp');
    colormap(ax_tA, 'jet'); clim(ax_tA, c_lim_topo);
    title(ax_tA, sprintf('%d ms', round(topo_times(k)*1000)), 'FontSize', 11, 'FontWeight', 'bold');
end

% Panel B Topoplots (150, 400, 900 ms)
for k = 1:length(topo_times)
    frac_x   = (topo_times(k) - time_win(1)) / t_span;
    center_x = ax2_pos(1) + frac_x * ax2_pos(3);
    ax_tB    = axes('Position', [center_x - topo_w/2, topo_y, topo_w, topo_h]); %#ok<LAXES>
    topoplot(topo_vecs_B(:, k), EEG.chanlocs, 'maplimits', c_lim_topo, ...
        'numcontour', 3, 'electrodes', 'on', 'style', 'map', 'shading', 'interp');
    colormap(ax_tB, 'jet'); clim(ax_tB, c_lim_topo);
    title(ax_tB, sprintf('%d ms', round(topo_times(k)*1000)), 'FontSize', 11, 'FontWeight', 'bold');
end

% Colorbar shifted left to x = 0.468 to eliminate overlap with ax2 y-axis label
ax_cb = axes('Position', [0.468, 0.48, 0.012, 0.18], 'Visible', 'off');
colormap(ax_cb, 'jet'); clim(ax_cb, c_lim_topo);
cb = colorbar(ax_cb, 'Position', [0.468, 0.48, 0.012, 0.18]);
cb.FontSize = 11;
cb.Title.String = cb_unit_str;
cb.Title.FontSize = 12;
cb.Title.FontWeight = 'bold';

% Panel Titles
ax_title = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
text(ax_title, ax1_pos(1) + ax1_pos(3)/2, 0.965, 'Stimulus Delivered (C3)', ...
    'HorizontalAlignment', 'center', 'FontSize', 22, 'FontWeight', 'bold');
text(ax_title, ax2_pos(1) + ax2_pos(3)/2, 0.965, 'Stimulus Omitted / Delayed (C3)', ...
    'HorizontalAlignment', 'center', 'FontSize', 22, 'FontWeight', 'bold');

% -------------------------------------------------------------------------
% 8. EXPORT 600 PPI PUBLICATION IMAGE
% -------------------------------------------------------------------------
fig1B_out_png = fullfile(output_path, 'Figure_1B_Channel_C3_Dynamics_600ppi.png');
fig1B_out_svg = fullfile(output_path, 'Figure_1B_Channel_C3_Dynamics.svg');

try
    exportgraphics(fig1B, fig1B_out_png, 'Resolution', 600);
    saveas(fig1B, fig1B_out_svg);
    fprintf('  -> Saved Figure 1B Visual: %s\n', fig1B_out_png);
    fprintf('  -> Saved Figure 1B Vector SVG: %s\n', fig1B_out_svg);
catch
    print(fig1B, fig1B_out_png, '-dpng', '-r600');
end
close(fig1B);
warning(orig_warn_1B);

%% =========================================================================
%% LOCAL HELPER FUNCTIONS FOR FIGURE 1B WILCOXON BARS
%% =========================================================================
function sig_mask = apply_bh_fdr(p_vals, alpha_lvl, do_fdr)
    p_vals = p_vals(:)';
    if ~do_fdr
        sig_mask = (p_vals < alpha_lvl);
        return;
    end
    valid = ~isnan(p_vals);
    p_sub = p_vals(valid);
    m     = length(p_sub);
    [p_sorted, ~] = sort(p_sub);
    thresh_vec    = (1:m) * (alpha_lvl / m);
     below_idx     = find(p_sorted <= thresh_vec, 1, 'last');
    if isempty(below_idx)
        sig_mask = false(size(p_vals));
    else
        p_crit   = p_sorted(below_idx);
        sig_mask = (p_vals <= p_crit);
    end
end

function clean_mask = filter_contiguous_runs(sig_mask, min_pts)
    clean_mask = false(size(sig_mask));
    d          = diff([0, double(sig_mask(:)'), 0]);
    run_starts = find(d == 1);
    run_ends   = find(d == -1) - 1;
    for r = 1:length(run_starts)
        if (run_ends(r) - run_starts(r) + 1) >= min_pts
            clean_mask(run_starts(r):run_ends(r)) = true;
        end
    end
end

function draw_sig_bars(ax_h, t_vec, sig_mask, y_level, bar_col, lw)
    d          = diff([0, double(sig_mask(:)'), 0]);
    run_starts = find(d == 1);
    run_ends   = find(d == -1) - 1;
    for r = 1:length(run_starts)
        t1 = t_vec(run_starts(r));
        t2 = t_vec(run_ends(r));
        plot(ax_h, [t1, t2], [y_level, y_level], '-', ...
            'Color', bar_col, 'LineWidth', lw, 'HandleVisibility', 'off');
    end
end

%% =========================================================================
%% FIGURE 1C: MULTI-CHANNEL STACKED ERPs (1x5 LAYOUT, [-0.6, 1.0 s] WINDOW)
%% =========================================================================
fprintf('\nGenerating Figure 1C: Multi-Channel Stacked ERPs ([-0.6, 1.0 s], Trial Average +/- 1 SD)...\n');

subj_to_plot = 1; 

plot_chan_labels = {'AF3', 'F3', 'FZ', 'FC5', 'FC1', 'T7', 'C3', 'Cz', ...
                    'C4', 'T8', 'TP7', 'CP3', 'CP1', 'P5', 'P1', 'O1', 'Oz', 'O2'};
n_plot_chans = length(plot_chan_labels);

time_win_sec = [-0.600, 1.000];

chan_indices = zeros(n_plot_chans, 1);
if exist('EEG', 'var') && isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs)
    eeg_labels = upper({EEG.chanlocs.labels});
    for k = 1:n_plot_chans
        idx = find(strcmpi(eeg_labels, upper(plot_chan_labels{k})), 1);
        if ~isempty(idx)
            chan_indices(k) = idx;
        else
            chan_indices(k) = min(k, cfg.num_ch);
        end
    end
else
    chan_indices = 1:min(n_plot_chans, cfg.num_ch);
end

time_mask = (time_s >= time_win_sec(1)) & (time_s <= time_win_sec(2));
t_plot    = time_s(time_mask);

iLocalBase = (t_plot >= -0.600) & (t_plot <= -0.500);

c_P1      = [0.07, 0.35, 0.37];
c_P2_500  = [0.47, 0.12, 0.08];
c_P3_500  = [0.88, 0.57, 0.10];
c_P2_2000 = [0.90, 0.59, 0.51];
c_P3_null = [0.94, 0.75, 0.47];

all_conds_1c = {'P1', 'P2_500', 'P3_500', 'P2_2000', 'P3_missing'};
names_1c     = {'Cued (P1)', 'Unpred 500 (P2-500)', 'Rand 500 (P3-500)', 'Unpred 2000 (P2-2000)', 'Rand Null (P3-null)'};
colors_1c    = {c_P1, c_P2_500, c_P3_500, c_P2_2000, c_P3_null};

data_1c  = struct();
all_stds = [];
all_p2p  = [];

for col = 1:length(all_conds_1c)
    c_name = all_conds_1c{col};
    
    if isnumeric(subj_to_plot)
        raw_mat   = data_all_conds.Raw.(c_name){subj_to_plot};
        sub_chans = raw_mat(chan_indices, time_mask, :);
        m_stack   = mean(sub_chans, 3, 'omitnan');
        s_stack   = std(sub_chans, 0, 3, 'omitnan');
    else
        subjs_valid = valid_retained_subjs;
        subj_means  = zeros(n_plot_chans, length(t_plot), length(subjs_valid));
        for s_idx = 1:length(subjs_valid)
            s_num = subjs_valid(s_idx);
            raw_s = data_all_conds.Raw.(c_name){s_num};
            subj_means(:, :, s_idx) = mean(raw_s(chan_indices, time_mask, :), 3, 'omitnan');
        end
        m_stack = mean(subj_means, 3, 'omitnan');
        s_stack = std(subj_means, 0, 3, 'omitnan');
    end
    
    if any(iLocalBase)
        m_stack = m_stack - mean(m_stack(:, iLocalBase), 2, 'omitnan');
    end
    
    data_1c.(c_name).mean = m_stack;
    data_1c.(c_name).std  = s_stack;
    
    all_stds = [all_stds; s_stack(:)];
    all_p2p  = [all_p2p; max(m_stack, [], 2) - min(m_stack, [], 2)];
end

med_sd  = median(all_stds(~isnan(all_stds)));
p90_erp = prctile(all_p2p(~isnan(all_p2p)), 90);
if isempty(med_sd) || isnan(med_sd) || med_sd <= 0, med_sd = 8; end
if isempty(p90_erp) || isnan(p90_erp) || p90_erp <= 0, p90_erp = 8; end

chan_step = max(1.35 * med_sd, 1.40 * p90_erp); 
y_offsets = (n_plot_chans:-1:1)' * chan_step;

fig1C = figure('Position', [50, 50, 2200, 1500], 'Color', 'w', ...
    'Name', 'Figure 1C: Multi-Channel Stacked ERPs');
t_layout = tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'normal');

for col = 1:5
    c_name  = all_conds_1c{col};
    c_color = colors_1c{col};
    c_title = names_1c{col};
    
    ax = nexttile(col); hold on; set(gca, 'FontSize', 22, 'LineWidth', 1.5);
    box off;
    
    xline(-0.5, 'k--', 'LineWidth', 1.4, 'Alpha', 0.75, 'HandleVisibility', 'off');
    
    if ismember(c_name, {'P1', 'P2_500', 'P3_500'})
        xline(0.0, 'Color', [0.15 0.65 0.25], 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    else
        xline(0.0, 'Color', [0.45 0.45 0.45], 'LineStyle', ':', 'LineWidth', 1.4, 'HandleVisibility', 'off');
    end
    
    m_mat = data_1c.(c_name).mean;
    s_mat = data_1c.(c_name).std;
    
    t_fill = [t_plot, fliplr(t_plot)];
    for ch = 1:n_plot_chans
        base_y = y_offsets(ch);
        y_m    = m_mat(ch, :) + base_y;
        y_s    = min(s_mat(ch, :), 1.15 * chan_step);
        
        fill(t_fill, [y_m + y_s, fliplr(y_m - y_s)], c_color, ...
            'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    for ch = 1:n_plot_chans
        base_y = y_offsets(ch);
        y_m    = m_mat(ch, :) + base_y;
        plot(t_plot, y_m, 'Color', c_color, 'LineWidth', 1.9);
    end
    
    xlim([-0.6, 1.0]);
    xticks([-0.5, 0, 0.5, 1.0]);
    ylim([0, (n_plot_chans + 1) * chan_step]);
    
    yticks(sort(y_offsets));
    if col == 1
        yticklabels(fliplr(plot_chan_labels));
        ylabel('Channels (stacked)', 'FontSize', 26, 'FontWeight', 'bold');
    else
        yticklabels({});
    end
    
    title(c_title, 'FontSize', 18, 'FontWeight', 'bold', 'Color', c_color);
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold');
    grid on; set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
end

if isnumeric(subj_to_plot)
    title(t_layout, sprintf('Subject %d: Multi-Channel ERP Dynamics (Trial Average \\pm 1 SD)', subj_to_plot), ...
        'FontSize', 28, 'FontWeight', 'bold');
    out_tag = sprintf('Subj%02d', subj_to_plot);
else
    title(t_layout, 'Grand Average: Multi-Channel ERP Dynamics (Mean \\pm 1 SD across Subjects)', ...
        'FontSize', 28, 'FontWeight', 'bold');
    out_tag = 'GrandAvg';
end

% -------------------------------------------------------------------------
% 6. EXPORT 600 PPI PUBLICATION IMAGES (PNG + VECTOR SVG)
% -------------------------------------------------------------------------
% Changed destination from screening_fig_dir to output_path
fig1C_out_png = fullfile(output_path, sprintf('Figure_1C_Stacked_Channels_%s_600ppi.png', out_tag));
fig1C_out_svg = fullfile(output_path, sprintf('Figure_1C_Stacked_Channels_%s.svg', out_tag));

try
    exportgraphics(fig1C, fig1C_out_png, 'Resolution', 600);
    saveas(fig1C, fig1C_out_svg);
    fprintf('  -> Saved Figure 1C High-Res PNG: %s\n', fig1C_out_png);
    fprintf('  -> Saved Figure 1C Vector SVG: %s\n', fig1C_out_svg);
catch
    print(fig1C, fig1C_out_png, '-dpng', '-r600');
end
close(fig1C);

%% =========================================================================
%% FIGURE 2: 4-ROW TOPOPLOTS (RAW, ALPHA, BETA, BETA/ALPHA) OVER [-0.2, 1.0 s]
%% =========================================================================
fprintf('\nGenerating Figure 2: 4-Row Spectral Power & Ratio Index Topoplots [-0.2 to 1.0 s]...\n');

step_fig2    = 0.20; 
win_starts_2 = -0.200 : step_fig2 : (1.000 - step_fig2);
t_centers_2  = win_starts_2 + step_fig2 / 2;
n_win_2      = length(win_starts_2);

subjs_fig2   = valid_retained_subjs;
n_subjs_2    = length(subjs_fig2);

conds_fig2   = {'P1', 'P3_missing', 'P2_2000'};
pow_groups   = struct();

for c_i = 1:length(conds_fig2)
    c_name = conds_fig2{c_i};
    pow_cell = cell(n_subjs_2, 3);
    
    for s_i = 1:n_subjs_2
        s = subjs_fig2(s_i);
        
        X_raw = data_all_conds.Raw.(c_name){s};
        X_a   = data_all_conds.Alpha.(c_name){s};
        X_b   = data_all_conds.Beta.(c_name){s};
        [nc, nt, ntr] = size(X_a);
        
        valid_tr  = find(~isnan(squeeze(X_a(1, 1, :))));
        P_raw_win = nan(nc, n_win_2, length(valid_tr));
        P_a_win   = nan(nc, n_win_2, length(valid_tr));
        P_b_win   = nan(nc, n_win_2, length(valid_tr));
        
        for tr_i = 1:length(valid_tr)
            tr = valid_tr(tr_i);
            
            env_raw_sq = abs(hilbert(double(X_raw(:, :, tr))')).^2;
            env_a_sq   = abs(hilbert(double(X_a(:, :, tr))')).^2;
            env_b_sq   = abs(hilbert(double(X_b(:, :, tr))')).^2;
            
            for w = 1:n_win_2
                w_mask = (time_s >= win_starts_2(w)) & (time_s < win_starts_2(w) + step_fig2);
                P_raw_win(:, w, tr_i) = mean(env_raw_sq(w_mask, :), 1, 'omitnan')';
                P_a_win(:, w, tr_i)   = mean(env_a_sq(w_mask, :), 1, 'omitnan')';
                P_b_win(:, w, tr_i)   = mean(env_b_sq(w_mask, :), 1, 'omitnan')';
            end
        end
        
        pow_cell{s_i, 1} = P_raw_win;
        pow_cell{s_i, 2} = P_a_win;
        pow_cell{s_i, 3} = P_b_win;
    end
    pow_groups.(c_name) = pow_cell;
end

% 1. Figure 2A: (Cued - Rand Null) / Cued (4 Rows)
plot_group_power_index_v2(pow_groups.P1, pow_groups.P3_missing, t_centers_2, ...
    'Cued', 'Rand Null', 'Cued vs Rand Null (-0.2 to 1.0s)', ...
    EEG.chanlocs, 1000, 0.05, output_path);

% 2. Figure 2B: (Cued - Unpred 2000) / Cued (4 Rows)
plot_group_power_index_v2(pow_groups.P1, pow_groups.P2_2000, t_centers_2, ...
    'Cued', 'Unpred 2000', 'Cued vs Unpred 2000 (-0.2 to 1.0s)', ...
    EEG.chanlocs, 1000, 0.05, output_path);

%% =========================================================================
%% AGGREGATE FIGURE 3B: MULTI-BAND SCREE PLOT
%% =========================================================================
disp('Generating Figure 3B: Aggregate Cumulative Variance Scree Plot...');

fig3B = figure('Position', [100, 100, 1200, 800], 'Name', 'Figure 3B: Scree Plot');
hold on; 
set(gca, 'FontSize', 22);

bands = {'Raw', 'Alpha', 'Beta', 'BetaAlphaRatio'};
band_labels = {'Raw Broadband', 'Alpha (8-12 Hz)', 'Beta (15-30 Hz)', 'Beta/Alpha Ratio'};
colors = {[0 0 0], [0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4940 0.1840 0.5560]};
markers = {'o', 's', '^', 'd'};

for b = 1:length(bands)
    if isfield(aggregate_spectra, bands{b})
        plot(1:cfg.num_ch, aggregate_spectra.(bands{b}), ['-' markers{b}], ...
             'Color', colors{b}, 'LineWidth', 2.5, 'MarkerSize', 8, ...
             'MarkerFaceColor', colors{b}, 'DisplayName', band_labels{b});
    end
end

xline(cfg.m, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
yline(85, 'r:', 'LineWidth', 2, 'DisplayName', '85% Variance Threshold');
plot(cfg.m, aggregate_spectra.Raw(cfg.m), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
text(cfg.m + 0.2, 70, sprintf('Truncation at m=%d\n(p < 0.05 consistency)', cfg.m), ...
    'FontSize', 22, 'Color', 'k', 'FontWeight', 'bold', 'BackgroundColor', 'w');

xline(2, ':', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2, 'HandleVisibility', 'off');
text(2.2, 50, 'Beta/Alpha intrinsic 2D manifold (p < 0.05)', ...
    'FontSize', 26, 'Color', [0.4940 0.1840 0.5560], 'FontWeight', 'bold', 'BackgroundColor', 'w');
    
text(cfg.m + 0.2, 70, sprintf('Truncation at m=%d\n(>85%% variance)', cfg.m), ...
    'FontSize', 26, 'Color', 'k', 'FontWeight', 'bold', 'BackgroundColor', 'w');
    
xlim([0.5, 15]); xticks(1:15); ylim([0, 100]); yticks(0:10:100);
xlabel('Number of Components (Spatial Directions)', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Cumulative Variance Explained (%)', 'FontSize', 24, 'FontWeight', 'bold');
title('Intrinsic Dimensionality Across Frequency Bands', 'FontSize', 28, 'FontWeight', 'bold');

lgd = legend('Location', 'southeast'); lgd.FontSize = 20; grid on;

save_file_3b = fullfile(output_path, 'Figure_3B_Aggregate_ScreePlot');
try pause(0.5); saveas(fig3B, [save_file_3b, '.svg']); saveas(fig3B, [save_file_3b, '.png']); catch, end
close(fig3B);

%% =========================================================================
%% AGGREGATE FIGURE 5A: SPLIT INTO (1) RAW (1x2) & (2) ALPHA/BETA (2x2)
%% =========================================================================
disp('Generating Figure 5A: Split Raw (1x2) and Alpha & Beta (2x2) log(r) Subspace Redistribution...');

plot_bands_filt = {'Alpha', 'Beta'};
band_colors_filt = containers.Map({'Alpha', 'Beta'}, ...
    {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]});
raw_color = [0.15 0.15 0.15];
cond_markers_A = {'o', 's', 'd'}; cond_styles_A  = {'-', '--', ':'};
cond_markers_B = {'^', 'v', 'd', 'p'}; cond_styles_B  = {'-', '--', ':', '-.'};

all_y_5a = [];
for g_chk = 1:length(subj_group_names)
    for b_chk = {'Raw', 'Alpha', 'Beta'}
        b_name = b_chk{1};
        if isfield(aggregate_log_r, b_name) && length(aggregate_log_r.(b_name).cohort) >= g_chk
            m_val  = aggregate_log_r.(b_name).cohort(g_chk).mean;
            se_val = aggregate_log_r.(b_name).cohort(g_chk).se;
            all_y_5a = [all_y_5a; m_val(:) + se_val(:); m_val(:) - se_val(:)];
        end
    end
end
if isempty(all_y_5a), all_y_5a = [-1.5, 1.5]; end
global_5a_ymin = min(-1.5, floor(min(all_y_5a) * 1.15 * 2) / 2);
global_5a_ymax = max( 1.5, ceil(max(all_y_5a)  * 1.15 * 2) / 2);
star_offset = (global_5a_ymax - global_5a_ymin) * 0.03;

for g = 1:length(subj_group_names)
    cohort_str = strrep(subj_group_names{g}, '_', ' ');
    
    % =====================================================================
    % FIGURE 5A (PART 1): RAW BROADBAND SIGNAL (1x2 GRID)
    % =====================================================================
    fig5A_Raw = figure('Position', [100, 100, 1300, 520], ...
        'Name', sprintf('Figure 5A (Raw): %s', cohort_str), 'Color', 'w');
    tiledlayout(1, 2, 'TileSpacing', 'normal', 'Padding', 'normal');
    
    % --- COL 1: STIMULUS DELIVERED - RAW BROADBAND ---
    nexttile; hold on; set(gca, 'FontSize', 16);
    yline(0, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    if isfield(aggregate_log_r, 'Raw') && length(aggregate_log_r.Raw.cohort) >= g
        grp_m  = aggregate_log_r.Raw.cohort(g).mean;
        grp_se = aggregate_log_r.Raw.cohort(g).se;
        for i = 1:length(cfg.group_A)
            c_idx = find(strcmp(cfg.all_conds, cfg.group_A{i}));
            if isempty(c_idx), continue; end
            x_jitter = (1:cfg.m) + (i - (length(cfg.group_A)+1)/2) * 0.08;
            
            errorbar(x_jitter, grp_m(c_idx, :), grp_se(c_idx, :), grp_se(c_idx, :), ...
                'LineStyle', cond_styles_A{i}, 'Marker', cond_markers_A{i}, ...
                'Color', raw_color, 'LineWidth', 2.5, 'MarkerSize', 8, ...
                'MarkerFaceColor', raw_color, 'DisplayName', sprintf('%s (Raw)', clean_name(cfg.group_A{i})));
            
            if isfield(aggregate_log_r.Raw.cohort(g), 'subj_data')
                s_data = aggregate_log_r.Raw.cohort(g).subj_data;
                for k = 1:cfg.m
                    dist = squeeze(s_data(c_idx, k, :));
                    if length(dist) > 1 && ~all(isnan(dist))
                        [~, p_val] = ttest(dist, 0);
                        if p_val < 0.05
                            star = '*'; if p_val < 0.01, star = '**'; end; if p_val < 0.001, star = '***'; end
                            y_pos = grp_m(c_idx, k) + sign(grp_m(c_idx, k) + eps) * (grp_se(c_idx, k) + star_offset);
                            v_align = 'bottom'; if grp_m(c_idx, k) < 0, v_align = 'top'; end
                            text(x_jitter(k), y_pos, star, 'Color', raw_color, 'FontSize', 18, 'FontWeight', 'bold', ...
                                'HorizontalAlignment', 'center', 'VerticalAlignment', v_align);
                        end
                    end
                end
            end
        end
    end
    grid on; xlim([0.5, cfg.m + 0.5]); xticks(1:cfg.m);
    xlabel('Spatial Direction Index', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('log(r_i) \pm SEM', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Stimulus Delivered — Raw\n(%s)', cohort_str), 'FontSize', 16, 'FontWeight', 'bold');
    lgd = legend('Location', 'best'); lgd.FontSize = 14;
    ylim([global_5a_ymin, global_5a_ymax]);
    
    % --- COL 2: STIMULUS OMITTED / DELAYED - RAW BROADBAND ---
    nexttile; hold on; set(gca, 'FontSize', 16);
    yline(0, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    if isfield(aggregate_log_r, 'Raw') && length(aggregate_log_r.Raw.cohort) >= g
        grp_m  = aggregate_log_r.Raw.cohort(g).mean;
        grp_se = aggregate_log_r.Raw.cohort(g).se;
        for i = 1:length(cfg.group_B)
            c_idx = find(strcmp(cfg.all_conds, cfg.group_B{i}));
            if isempty(c_idx), continue; end
            x_jitter = (1:cfg.m) + (i - (length(cfg.group_B)+1)/2) * 0.08;
            
            errorbar(x_jitter, grp_m(c_idx, :), grp_se(c_idx, :), grp_se(c_idx, :), ...
                'LineStyle', cond_styles_B{i}, 'Marker', cond_markers_B{i}, ...
                'Color', raw_color, 'LineWidth', 2.5, 'MarkerSize', 8, ...
                'MarkerFaceColor', raw_color, 'DisplayName', sprintf('%s (Raw)', clean_name(cfg.group_B{i})));
            
            if isfield(aggregate_log_r.Raw.cohort(g), 'subj_data')
                s_data = aggregate_log_r.Raw.cohort(g).subj_data;
                for k = 1:cfg.m
                    dist = squeeze(s_data(c_idx, k, :));
                    if length(dist) > 1 && ~all(isnan(dist))
                        [~, p_val] = ttest(dist, 0);
                        if p_val < 0.05
                            star = '*'; if p_val < 0.01, star = '**'; end; if p_val < 0.001, star = '***'; end
                            y_pos = grp_m(c_idx, k) + sign(grp_m(c_idx, k) + eps) * (grp_se(c_idx, k) + star_offset);
                            v_align = 'bottom'; if grp_m(c_idx, k) < 0, v_align = 'top'; end
                            text(x_jitter(k), y_pos, star, 'Color', raw_color, 'FontSize', 18, 'FontWeight', 'bold', ...
                                'HorizontalAlignment', 'center', 'VerticalAlignment', v_align);
                        end
                    end
                end
            end
        end
    end
    grid on; xlim([0.5, cfg.m + 0.5]); xticks(1:cfg.m);
    xlabel('Spatial Direction Index', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Stimulus Omitted / Delayed — Raw\n(%s)', cohort_str), 'FontSize', 16, 'FontWeight', 'bold');
    lgd = legend('Location', 'best'); lgd.FontSize = 14;
    ylim([global_5a_ymin, global_5a_ymax]);
    
    sgtitle(sprintf('Subspace Redistribution Geometry — Raw (%s)', cohort_str), 'FontSize', 20, 'FontWeight', 'bold');
    
    save_file_5a_raw = fullfile(output_path, sprintf('Figure_5A_Raw_LogR_%s', fig_tags{g}));
    try pause(0.5); saveas(fig5A_Raw, [save_file_5a_raw, '.svg']); saveas(fig5A_Raw, [save_file_5a_raw, '.png']); catch, end
    close(fig5A_Raw);
    
    % =====================================================================
    % FIGURE 5A (PART 2): FILTERED ALPHA & BETA SIGNALS (2x2 GRID)
    % =====================================================================
    fig5A_AB = figure('Position', [100, 100, 1300, 920], ...
        'Name', sprintf('Figure 5A (Alpha & Beta): %s', cohort_str), 'Color', 'w');
    tiledlayout(2, 2, 'TileSpacing', 'normal', 'Padding', 'normal');
    
    % --- ROW 1, COL 1: STIMULUS DELIVERED - ALPHA ---
    b_name = 'Alpha'; col = band_colors_filt(b_name);
    nexttile; hold on; set(gca, 'FontSize', 16);
    yline(0, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    if isfield(aggregate_log_r, b_name) && length(aggregate_log_r.(b_name).cohort) >= g
        grp_m  = aggregate_log_r.(b_name).cohort(g).mean;
        grp_se = aggregate_log_r.(b_name).cohort(g).se;
        for i = 1:length(cfg.group_A)
            c_idx = find(strcmp(cfg.all_conds, cfg.group_A{i}));
            if isempty(c_idx), continue; end
            x_jitter = (1:cfg.m) + (i - (length(cfg.group_A)+1)/2) * 0.08;
            
            errorbar(x_jitter, grp_m(c_idx, :), grp_se(c_idx, :), grp_se(c_idx, :), ...
                'LineStyle', cond_styles_A{i}, 'Marker', cond_markers_A{i}, ...
                'Color', col, 'LineWidth', 2.5, 'MarkerSize', 8, ...
                'MarkerFaceColor', col, 'DisplayName', sprintf('%s (%s)', clean_name(cfg.group_A{i}), b_name));
            
            if isfield(aggregate_log_r.(b_name).cohort(g), 'subj_data')
                s_data = aggregate_log_r.(b_name).cohort(g).subj_data;
                for k = 1:cfg.m
                    dist = squeeze(s_data(c_idx, k, :));
                    if length(dist) > 1 && ~all(isnan(dist))
                        [~, p_val] = ttest(dist, 0);
                        if p_val < 0.05
                            star = '*'; if p_val < 0.01, star = '**'; end; if p_val < 0.001, star = '***'; end
                            y_pos = grp_m(c_idx, k) + sign(grp_m(c_idx, k) + eps) * (grp_se(c_idx, k) + star_offset);
                            v_align = 'bottom'; if grp_m(c_idx, k) < 0, v_align = 'top'; end
                            text(x_jitter(k), y_pos, star, 'Color', col, 'FontSize', 18, 'FontWeight', 'bold', ...
                                'HorizontalAlignment', 'center', 'VerticalAlignment', v_align);
                        end
                    end
                end
            end
        end
    end
    grid on; xlim([0.5, cfg.m + 0.5]); xticks(1:cfg.m); xticklabels({});
    ylabel('log(r_i) \pm SEM', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Stimulus Delivered — Alpha\n(%s)', cohort_str), 'FontSize', 16, 'FontWeight', 'bold');
    lgd = legend('Location', 'best'); lgd.FontSize = 14;
    ylim([global_5a_ymin, global_5a_ymax]);
    
    % --- ROW 1, COL 2: STIMULUS OMITTED / DELAYED - ALPHA ---
    nexttile; hold on; set(gca, 'FontSize', 16);
    yline(0, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    if isfield(aggregate_log_r, b_name) && length(aggregate_log_r.(b_name).cohort) >= g
        grp_m  = aggregate_log_r.(b_name).cohort(g).mean;
        grp_se = aggregate_log_r.(b_name).cohort(g).se;
        for i = 1:length(cfg.group_B)
            c_idx = find(strcmp(cfg.all_conds, cfg.group_B{i}));
            if isempty(c_idx), continue; end
            x_jitter = (1:cfg.m) + (i - (length(cfg.group_B)+1)/2) * 0.08;
            
            errorbar(x_jitter, grp_m(c_idx, :), grp_se(c_idx, :), grp_se(c_idx, :), ...
                'LineStyle', cond_styles_B{i}, 'Marker', cond_markers_B{i}, ...
                'Color', col, 'LineWidth', 2.5, 'MarkerSize', 8, ...
                'MarkerFaceColor', col, 'DisplayName', sprintf('%s (%s)', clean_name(cfg.group_B{i}), b_name));
            
            if isfield(aggregate_log_r.(b_name).cohort(g), 'subj_data')
                s_data = aggregate_log_r.(b_name).cohort(g).subj_data;
                for k = 1:cfg.m
                    dist = squeeze(s_data(c_idx, k, :));
                    if length(dist) > 1 && ~all(isnan(dist))
                        [~, p_val] = ttest(dist, 0);
                        if p_val < 0.05
                            star = '*'; if p_val < 0.01, star = '**'; end; if p_val < 0.001, star = '***'; end
                            y_pos = grp_m(c_idx, k) + sign(grp_m(c_idx, k) + eps) * (grp_se(c_idx, k) + star_offset);
                            v_align = 'bottom'; if grp_m(c_idx, k) < 0, v_align = 'top'; end
                            text(x_jitter(k), y_pos, star, 'Color', col, 'FontSize', 18, 'FontWeight', 'bold', ...
                                'HorizontalAlignment', 'center', 'VerticalAlignment', v_align);
                        end
                    end
                end
            end
        end
    end
    grid on; xlim([0.5, cfg.m + 0.5]); xticks(1:cfg.m); xticklabels({});
    title(sprintf('Stimulus Omitted / Delayed — Alpha\n(%s)', cohort_str), 'FontSize', 16, 'FontWeight', 'bold');
    lgd = legend('Location', 'best'); lgd.FontSize = 14;
    ylim([global_5a_ymin, global_5a_ymax]);
    
    % --- ROW 2, COL 1: STIMULUS DELIVERED - BETA ---
    b_name = 'Beta'; col = band_colors_filt(b_name);
    nexttile; hold on; set(gca, 'FontSize', 16);
    yline(0, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    if isfield(aggregate_log_r, b_name) && length(aggregate_log_r.(b_name).cohort) >= g
        grp_m  = aggregate_log_r.(b_name).cohort(g).mean;
        grp_se = aggregate_log_r.(b_name).cohort(g).se;
        for i = 1:length(cfg.group_A)
            c_idx = find(strcmp(cfg.all_conds, cfg.group_A{i}));
            if isempty(c_idx), continue; end
            x_jitter = (1:cfg.m) + (i - (length(cfg.group_A)+1)/2) * 0.08;
            
            errorbar(x_jitter, grp_m(c_idx, :), grp_se(c_idx, :), grp_se(c_idx, :), ...
                'LineStyle', cond_styles_A{i}, 'Marker', cond_markers_A{i}, ...
                'Color', col, 'LineWidth', 2.5, 'MarkerSize', 8, ...
                'MarkerFaceColor', col, 'DisplayName', sprintf('%s (%s)', clean_name(cfg.group_A{i}), b_name));
            
            if isfield(aggregate_log_r.(b_name).cohort(g), 'subj_data')
                s_data = aggregate_log_r.(b_name).cohort(g).subj_data;
                for k = 1:cfg.m
                    dist = squeeze(s_data(c_idx, k, :));
                    if length(dist) > 1 && ~all(isnan(dist))
                        [~, p_val] = ttest(dist, 0);
                        if p_val < 0.05
                            star = '*'; if p_val < 0.01, star = '**'; end; if p_val < 0.001, star = '***'; end
                            y_pos = grp_m(c_idx, k) + sign(grp_m(c_idx, k) + eps) * (grp_se(c_idx, k) + star_offset);
                            v_align = 'bottom'; if grp_m(c_idx, k) < 0, v_align = 'top'; end
                            text(x_jitter(k), y_pos, star, 'Color', col, 'FontSize', 18, 'FontWeight', 'bold', ...
                                'HorizontalAlignment', 'center', 'VerticalAlignment', v_align);
                        end
                    end
                end
            end
        end
    end
    grid on; xlim([0.5, cfg.m + 0.5]); xticks(1:cfg.m);
    xlabel('Spatial Direction Index', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('log(r_i) \pm SEM', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Stimulus Delivered — Beta\n(%s)', cohort_str), 'FontSize', 16, 'FontWeight', 'bold');
    lgd = legend('Location', 'best'); lgd.FontSize = 14;
    ylim([global_5a_ymin, global_5a_ymax]);
    
    % --- ROW 2, COL 2: STIMULUS OMITTED / DELAYED - BETA ---
    nexttile; hold on; set(gca, 'FontSize', 16);
    yline(0, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    if isfield(aggregate_log_r, b_name) && length(aggregate_log_r.(b_name).cohort) >= g
        grp_m  = aggregate_log_r.(b_name).cohort(g).mean;
        grp_se = aggregate_log_r.(b_name).cohort(g).se;
        for i = 1:length(cfg.group_B)
            c_idx = find(strcmp(cfg.all_conds, cfg.group_B{i}));
            if isempty(c_idx), continue; end
            x_jitter = (1:cfg.m) + (i - (length(cfg.group_B)+1)/2) * 0.08;
            
            errorbar(x_jitter, grp_m(c_idx, :), grp_se(c_idx, :), grp_se(c_idx, :), ...
                'LineStyle', cond_styles_B{i}, 'Marker', cond_markers_B{i}, ...
                'Color', col, 'LineWidth', 2.5, 'MarkerSize', 8, ...
                'MarkerFaceColor', col, 'DisplayName', sprintf('%s (%s)', clean_name(cfg.group_B{i}), b_name));
            
            if isfield(aggregate_log_r.(b_name).cohort(g), 'subj_data')
                s_data = aggregate_log_r.(b_name).cohort(g).subj_data;
                for k = 1:cfg.m
                    dist = squeeze(s_data(c_idx, k, :));
                    if length(dist) > 1 && ~all(isnan(dist))
                        [~, p_val] = ttest(dist, 0);
                        if p_val < 0.05
                            star = '*'; if p_val < 0.01, star = '**'; end; if p_val < 0.001, star = '***'; end
                            y_pos = grp_m(c_idx, k) + sign(grp_m(c_idx, k) + eps) * (grp_se(c_idx, k) + star_offset);
                            v_align = 'bottom'; if grp_m(c_idx, k) < 0, v_align = 'top'; end
                            text(x_jitter(k), y_pos, star, 'Color', col, 'FontSize', 18, 'FontWeight', 'bold', ...
                                'HorizontalAlignment', 'center', 'VerticalAlignment', v_align);
                        end
                    end
                end
            end
        end
    end
    grid on; xlim([0.5, cfg.m + 0.5]); xticks(1:cfg.m);
    xlabel('Spatial Direction Index', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Stimulus Omitted / Delayed — Beta\n(%s)', cohort_str), 'FontSize', 16, 'FontWeight', 'bold');
    lgd = legend('Location', 'best'); lgd.FontSize = 14;
    ylim([global_5a_ymin, global_5a_ymax]);
    
    sgtitle(sprintf('Subspace Redistribution Geometry — Alpha & Beta (%s)', cohort_str), 'FontSize', 20, 'FontWeight', 'bold');
    
    save_file_5a_ab = fullfile(output_path, sprintf('Figure_5A_AlphaBeta_LogR_%s', fig_tags{g}));
    try pause(0.5); saveas(fig5A_AB, [save_file_5a_ab, '.svg']); saveas(fig5A_AB, [save_file_5a_ab, '.png']); catch, end
    close(fig5A_AB);
end

%% =========================================================================
%% AGGREGATE FIGURE 5B: SPLIT GLOBAL GAIN (G) BY STIMULUS DELIVERED & OMITTED / DELAYED
%% =========================================================================
disp('Generating Figure 5B: Split Stimulus Delivered & Stimulus Omitted / Delayed Gain Plots (Wilcoxon Signed-Rank)...');
plot_bands  = {'Raw', 'Alpha', 'Beta', 'BetaAlphaRatio'};
band_labels = {'Raw', 'Alpha', 'Beta', sprintf('Beta/Alpha\nRatio')};
c_map = containers.Map();
c_map('P1')         = [0.500 0.500 0.500];
c_map('P2_500')     = [0.850 0.325 0.098];
c_map('P3_500')     = [0.929 0.694 0.125];
c_map('P2_2000')    = [0.466 0.674 0.188];
c_map('P3_missing') = [0.301 0.745 0.933];
conds_plot_A = [{'P1'}, cfg.group_A];                  
conds_plot_B = [{'P1'}, intersect(cfg.group_B, cfg.all_conds, 'stable')]; 
idx_P1 = find(strcmp(cfg.all_conds, 'P1'));
global_max_gain = 0;
for g_chk = 1:length(subj_group_names)
    for b = 1:length(plot_bands)
        b_name = plot_bands{b};
        if isfield(aggregate_G, b_name) && length(aggregate_G.(b_name).cohort) >= g_chk
            m_v  = aggregate_G.(b_name).cohort(g_chk).mean;
            se_v = aggregate_G.(b_name).cohort(g_chk).se;
            global_max_gain = max(global_max_gain, max(m_v + se_v));
        end
    end
end
global_y_lim_5b = max(global_max_gain * 1.15, 2.5);

for g = 1:length(subj_group_names)
    cohort_str = strrep(subj_group_names{g}, '_', ' ');
    fig5B = figure('Position', [100, 100, 1750, 800], 'Name', sprintf('Figure 5B: Gain %s', cohort_str));
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    % SUBPLOT 1: STIMULUS DELIVERED
    nexttile; hold on; set(gca, 'FontSize', 18);
    
    nA = length(conds_plot_A);
    Y_mean_A = zeros(length(plot_bands), nA);
    Y_se_A   = zeros(length(plot_bands), nA);
    
    for b = 1:length(plot_bands)
        b_name = plot_bands{b};
        if isfield(aggregate_G, b_name) && length(aggregate_G.(b_name).cohort) >= g
            for i = 1:nA
                c_idx = find(strcmp(cfg.all_conds, conds_plot_A{i}));
                if ~isempty(c_idx)
                    Y_mean_A(b, i) = aggregate_G.(b_name).cohort(g).mean(c_idx);
                    Y_se_A(b, i)   = aggregate_G.(b_name).cohort(g).se(c_idx);
                end
            end
        end
    end
    
    b_hA = bar(Y_mean_A, 'grouped');
    for i = 1:nA
        b_hA(i).FaceColor = c_map(conds_plot_A{i});
        b_hA(i).DisplayName = clean_name(conds_plot_A{i});
    end
    
    groupwidthA = min(0.8, nA / (nA + 1.5));
    get_xA = @(b, i) b - groupwidthA/2 + (2*i - 1) * groupwidthA / (2*nA);
    for i = 1:nA
        x = (1:length(plot_bands)) - groupwidthA/2 + (2*i - 1) * groupwidthA / (2*nA);
        errorbar(x, Y_mean_A(:, i), Y_se_A(:, i), 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    yline(1.0, 'k--', 'LineWidth', 2, 'DisplayName', 'Cued Baseline (1.0)');
    
    for b = 1:length(plot_bands)
        b_name = plot_bands{b};
        if ~isfield(aggregate_G, b_name) || length(aggregate_G.(b_name).cohort) < g, continue; end
        G_data = aggregate_G.(b_name).cohort(g).subj_data;
        if size(G_data, 2) <= 1, continue; end
        
        band_max = max(Y_mean_A(b, :) + Y_se_A(b, :));
        y_curr = band_max + global_y_lim_5b * 0.05;
        
        for i = 2:nA
            c_idx = find(strcmp(cfg.all_conds, conds_plot_A{i}));
            if isempty(c_idx) || isempty(idx_P1), continue; end
            
            p_right = signrank(G_data(c_idx, :), G_data(idx_P1, :), 'tail', 'right');
            p_left  = signrank(G_data(c_idx, :), G_data(idx_P1, :), 'tail', 'left');
            p_val = min(p_right, p_left);
            
            if p_val < 0.05
                x1 = get_xA(b, 1); x2 = get_xA(b, i);
                h_bar = global_y_lim_5b * 0.02;
                plot([x1, x1, x2, x2], [y_curr, y_curr + h_bar, y_curr + h_bar, y_curr], 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
                star = '*'; if p_val < 0.01, star = '**'; end; if p_val < 0.001, star = '***'; end
                text(mean([x1, x2]), y_curr + h_bar * 1.5, star, 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                y_curr = y_curr + global_y_lim_5b * 0.08;
            end
        end
    end
    
    xticks(1:length(plot_bands)); xticklabels(band_labels);
    xtickangle(0);
    ylabel('Global Subspace Gain (G)', 'FontSize', 20, 'FontWeight', 'bold');
    title(sprintf('Stimulus Delivered Gain (%s)', cohort_str), 'FontSize', 20, 'FontWeight', 'bold');
    ylim([0, global_y_lim_5b]); grid on;
    lgd = legend('Location', 'north'); lgd.FontSize = 14;
    
    % SUBPLOT 2: STIMULUS OMITTED / DELAYED
    nexttile; hold on; set(gca, 'FontSize', 18);
    
    nB = length(conds_plot_B);
    Y_mean_B = zeros(length(plot_bands), nB);
    Y_se_B   = zeros(length(plot_bands), nB);
    
    for b = 1:length(plot_bands)
        b_name = plot_bands{b};
        if isfield(aggregate_G, b_name) && length(aggregate_G.(b_name).cohort) >= g
            for i = 1:nB
                c_idx = find(strcmp(cfg.all_conds, conds_plot_B{i}));
                if ~isempty(c_idx)
                    Y_mean_B(b, i) = aggregate_G.(b_name).cohort(g).mean(c_idx);
                    Y_se_B(b, i)   = aggregate_G.(b_name).cohort(g).se(c_idx);
                end
            end
        end
    end
    
    b_hB = bar(Y_mean_B, 'grouped');
    for i = 1:nB
        b_hB(i).FaceColor = c_map(conds_plot_B{i});
        b_hB(i).DisplayName = clean_name(conds_plot_B{i});
    end
    
    groupwidthB = min(0.8, nB / (nB + 1.5));
    get_xB = @(b, i) b - groupwidthB/2 + (2*i - 1) * groupwidthB / (2*nB);
    for i = 1:nB
        x = (1:length(plot_bands)) - groupwidthB/2 + (2*i - 1) * groupwidthB / (2*nB);
        errorbar(x, Y_mean_B(:, i), Y_se_B(:, i), 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    yline(1.0, 'k--', 'LineWidth', 2, 'DisplayName', 'Cued Baseline (1.0)');
    
    for b = 1:length(plot_bands)
        b_name = plot_bands{b};
        if ~isfield(aggregate_G, b_name) || length(aggregate_G.(b_name).cohort) < g, continue; end
        G_data = aggregate_G.(b_name).cohort(g).subj_data;
        if size(G_data, 2) <= 1, continue; end
        
        band_max = max(Y_mean_B(b, :) + Y_se_B(b, :));
        y_curr = band_max + global_y_lim_5b * 0.05;
        
        for i = 2:nB
            c_idx = find(strcmp(cfg.all_conds, conds_plot_B{i}));
            if isempty(c_idx) || isempty(idx_P1), continue; end
            
            p_right = signrank(G_data(c_idx, :), G_data(idx_P1, :), 'tail', 'right');
            p_left  = signrank(G_data(c_idx, :), G_data(idx_P1, :), 'tail', 'left');
            p_val = min(p_right, p_left);
            
            if p_val < 0.05
                x1 = get_xB(b, 1); x2 = get_xB(b, i);
                h_bar = global_y_lim_5b * 0.02;
                plot([x1, x1, x2, x2], [y_curr, y_curr + h_bar, y_curr + h_bar, y_curr], 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
                star = '*'; if p_val < 0.01, star = '**'; end; if p_val < 0.001, star = '***'; end
                text(mean([x1, x2]), y_curr + h_bar * 1.5, star, 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                y_curr = y_curr + global_y_lim_5b * 0.08;
            end
        end
    end
    
    xticks(1:length(plot_bands)); xticklabels(band_labels);
    xtickangle(0);
    title(sprintf('Stimulus Omitted / Delayed Gain (%s)', cohort_str), 'FontSize', 20, 'FontWeight', 'bold');
    ylim([0, global_y_lim_5b]); grid on;
    lgd = legend('Location', 'north'); lgd.FontSize = 14;
    
    sgtitle(sprintf('Global Subspace Gain Scaling (%s)', cohort_str), 'FontSize', 24, 'FontWeight', 'bold');
    
    save_file_5b = fullfile(output_path, sprintf('Figure_5B_MultiBand_Gain_%s', fig_tags{g}));
    try pause(0.5); saveas(fig5B, [save_file_5b, '.svg']); saveas(fig5B, [save_file_5b, '.png']); catch, end
    close(fig5B);
end

%% =========================================================================
%% AGGREGATE FIGURE 5C: MULTI-BAND DELTA LATENCY SCATTER PLOT
%% =========================================================================
disp('Generating Figure 5C: Multi-Band Delta Latency Scatter Plots...');
plot_bands  = {'Raw', 'Alpha', 'Beta', 'BetaAlphaRatio'};
band_titles = {'Raw Broadband', 'Alpha (8-12 Hz)', 'Beta (15-30 Hz)', 'Beta/Alpha Ratio'};
colors_conds = [0.500 0.500 0.500; ...
                0.850 0.325 0.098; ...
                0.929 0.694 0.125; ...
                0.466 0.674 0.188; ...
                0.301 0.745 0.933; ...
                0.494 0.184 0.556];
erp_conds = {'P1', 'P2_500', 'P3_500'};
erp_conds_idx = find(ismember(cfg.all_conds, erp_conds));
for g = 1:length(subj_group_names)
    cohort_str = strrep(subj_group_names{g}, '_', ' ');
    fig5C = figure('Position', [100, 100, 1600, 700], 'Name', sprintf('Figure 5C: Latency %s', cohort_str));
    tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    all_lat = [];
    for b = 1:length(plot_bands)
        if isfield(aggregate_lat, plot_bands{b}) && length(aggregate_lat.(plot_bands{b}).cohort) >= g
            lat_subset = aggregate_lat.(plot_bands{b}).cohort(g).mean(erp_conds_idx, :);
            all_lat = [all_lat; lat_subset(:)];
        end
    end
    
    if isempty(all_lat) || all(isnan(all_lat)), all_lat = [-50, 50]; end
    y_bound = max(50, ceil(max(abs(all_lat)) / 25) * 25); 
    for b = 1:length(plot_bands)
        b_name = plot_bands{b};
        if ~isfield(aggregate_lat, b_name) || length(aggregate_lat.(b_name).cohort) < g, continue; end
        
        lat_m = aggregate_lat.(b_name).cohort(g).mean;
        
        nexttile; hold on;
        set(gca, 'FontSize', 18);
        
        yline(0, 'k-', 'LineWidth', 2.5, 'HandleVisibility', 'off');
        
        for i = 1:length(erp_conds_idx)
            c_idx = erp_conds_idx(i);
            
            rng(42); 
            x_jitter = i + (rand(1, cfg.m) - 0.5) * 0.35; 
            
            scatter(x_jitter, lat_m(c_idx, :), 100, colors_conds(c_idx, :), 'filled', ...
                'MarkerEdgeColor', 'k', 'DisplayName', clean_name(cfg.all_conds{c_idx}));
                
            mean_val = mean(lat_m(c_idx, :), 'omitnan');
            plot([i-0.3, i+0.3], [mean_val, mean_val], '-', 'Color', colors_conds(c_idx, :), 'LineWidth', 4, 'HandleVisibility', 'off');
            
            if c_idx ~= find(strcmp(cfg.all_conds, 'P1'))
                % Subject-level latency shifts across the n=12 retained subjects
                lat_data = aggregate_lat.(b_name).cohort(g).subj_data(c_idx, :);
                lat_data_clean = lat_data(~isnan(lat_data));
                if length(lat_data_clean) > 1
                    p_val = signrank(lat_data_clean, 0); % Non-parametric Wilcoxon signed-rank test
                else
                    p_val = 1;
                end
                if p_val < 0.05
                    star = '*'; if p_val < 0.01, star = '**'; end; if p_val < 0.001, star = '***'; end
                    test_mean = mean(lat_m(c_idx, :), 'omitnan');
                    if test_mean > 0
                        y_star = max(lat_m(c_idx, :)) + y_bound * 0.08;
                    else
                        y_star = min(lat_m(c_idx, :)) - y_bound * 0.08;
                    end
                    text(i, y_star, star, 'FontSize', 24, 'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'center', 'Color', colors_conds(c_idx, :));
                end
            end
        end
        
        xlim([0.5, 3.5]);
        xticks(1:3);
        xticklabels({'Cued', 'Unpred 500', 'Rand 500'});
        xtickangle(45);
        
        if b == 1, ylabel('\Delta Latency relative to Cued (ms)', 'FontSize', 22, 'FontWeight', 'bold'); end
        title(band_titles{b}, 'FontSize', 22, 'FontWeight', 'bold');
        ylim([-y_bound * 1.15, y_bound * 1.15]); 
        grid on;
    end
    
    sgtitle(sprintf('Stimulus Delivered: Tactile ERP \\Delta Latency Shifts (%s)', cohort_str), 'FontSize', 28, 'FontWeight', 'bold');
    
    lgd_names = cellfun(clean_name, cfg.all_conds(erp_conds_idx), 'UniformOutput', false);
    lgd = legend(lgd_names, 'Location', 'northeastoutside');
    lgd.FontSize = 18;
    
    save_file_5c = fullfile(output_path, sprintf('Figure_5C_MultiBand_Latency_%s', fig_tags{g}));
    try pause(0.5); saveas(fig5C, [save_file_5c, '.svg']); saveas(fig5C, [save_file_5c, '.png']); catch, end
    close(fig5C);
end

%% =========================================================================
%% AGGREGATE FIGURE 5D: MULTI-BAND % SUSTAINED (TONIC VS PHASIC) BAR PLOTS
%% =========================================================================
disp('Generating Figure 5D: Multi-Band Sustained % Bar Plots (Wilcoxon Signed-Rank)...');
idx_P1      = find(strcmp(cfg.all_conds, 'P1'));
idx_P2_500  = find(strcmp(cfg.all_conds, 'P2_500'));
idx_P3_500  = find(strcmp(cfg.all_conds, 'P3_500'));
idx_P2_2000 = find(strcmp(cfg.all_conds, 'P2_2000'));
idx_P3_mis  = find(strcmp(cfg.all_conds, 'P3_missing'));
for g = 1:length(subj_group_names)
    cohort_str = strrep(subj_group_names{g}, '_', ' ');
    fig5D = figure('Position', [100, 100, 1400, 800], 'Name', sprintf('Figure 5D: Sustained %s', cohort_str));
    hold on; set(gca, 'FontSize', 22);
    
    Y_mean = zeros(length(plot_bands), length(cfg.all_conds));
    Y_se   = zeros(length(plot_bands), length(cfg.all_conds));
    
    for b = 1:length(plot_bands)
        b_name = plot_bands{b};
        if isfield(aggregate_sust, b_name) && length(aggregate_sust.(b_name).cohort) >= g
            Y_mean(b, :) = reshape(aggregate_sust.(b_name).cohort(g).mean, 1, []);
            Y_se(b, :)   = reshape(aggregate_sust.(b_name).cohort(g).se, 1, []);
        end
    end
    
    b_handle = bar(Y_mean, 'grouped');
    
    for i = 1:length(b_handle)
        b_handle(i).FaceColor = colors_conds(i, :);
        b_handle(i).DisplayName = clean_name(cfg.all_conds{i});
    end
    
    ngroups = size(Y_mean, 1);
    nbars = size(Y_mean, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    get_x_pos = @(b, c) b - groupwidth/2 + (2*c-1) * groupwidth / (2*nbars);
    
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, Y_mean(:,i), Y_se(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    
    xticks(1:length(plot_bands));
    xticklabels({'Raw', 'Alpha', 'Beta', 'Beta/Alpha Ratio'});
    ylabel('Sustained Variance (%)', 'FontSize', 24, 'FontWeight', 'bold');
    title(sprintf('Sustained vs. Phasic Stability (%s)', cohort_str), 'FontSize', 28, 'FontWeight', 'bold');
    
    lgd = legend('Location', 'north');
    lgd.FontSize = 18;
    grid on;
    
    ylim([0, 60]);
    
    plot([1.6, 3.4], [3.5, 3.5], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot([1.6, 1.6], [2.5, 3.5], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot([3.4, 3.4], [2.5, 3.5], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    text(2.5, 5.0, 'Purely Phasic (< 3%)', 'FontSize', 20, 'Color', 'k', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'Margin', 2);
    
    % 1. Raw Band (b=1): Cued vs Rand 500 & Cued vs Rand Missing
    if isfield(aggregate_sust, 'Raw') && length(aggregate_sust.Raw.cohort) >= g && ~isempty(idx_P1)
        raw_data = aggregate_sust.Raw.cohort(g).subj_data;
        if size(raw_data, 2) > 1 
            y_shift = 4; 
            
            if ~isempty(idx_P3_500)
                p_r1 = signrank(raw_data(idx_P1, :), raw_data(idx_P3_500, :), 'tail', 'right');
                p_l1 = signrank(raw_data(idx_P1, :), raw_data(idx_P3_500, :), 'tail', 'left');
                p1   = min(p_r1, p_l1);
                
                if ~isnan(p1) && p1 < 0.05
                    x1 = get_x_pos(1, idx_P1); x2 = get_x_pos(1, idx_P3_500);
                    y = max(Y_mean(1, [idx_P1, idx_P3_500]) + Y_se(1, [idx_P1, idx_P3_500])) + y_shift;
                    plot([x1, x1, x2, x2], [y, y+1.5, y+1.5, y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    star = '*'; if p1 < 0.01, star = '**'; end; if p1 < 0.001, star = '***'; end
                    text(mean([x1, x2]), y+3, star, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                    y_shift = y_shift + 8; 
                end
            end
            
            if ~isempty(idx_P3_mis)
                p_r2 = signrank(raw_data(idx_P1, :), raw_data(idx_P3_mis, :), 'tail', 'right');
                p_l2 = signrank(raw_data(idx_P1, :), raw_data(idx_P3_mis, :), 'tail', 'left');
                p2   = min(p_r2, p_l2);
                
                if ~isnan(p2) && p2 < 0.05
                    x1 = get_x_pos(1, idx_P1); x2 = get_x_pos(1, idx_P3_mis);
                    y = max(Y_mean(1, [idx_P1, idx_P3_mis]) + Y_se(1, [idx_P1, idx_P3_mis])) + y_shift;
                    plot([x1, x1, x2, x2], [y, y+1.5, y+1.5, y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    star = '*'; if p2 < 0.01, star = '**'; end; if p2 < 0.001, star = '***'; end
                    text(mean([x1, x2]), y+3, star, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                end
            end
        end
    end
    
    % 2. Alpha (b=2) and Beta (b=3) Bands: Cued vs Rand 500 & Cued vs Rand Missing
    for b_idx = 2:3
        b_name = plot_bands{b_idx};
        if isfield(aggregate_sust, b_name) && length(aggregate_sust.(b_name).cohort) >= g && ~isempty(idx_P1)
            sust_data = aggregate_sust.(b_name).cohort(g).subj_data;
            if size(sust_data, 2) > 1 
                y_shift = 8; 
                
                if ~isempty(idx_P3_500)
                    p_rab1 = signrank(sust_data(idx_P1, :), sust_data(idx_P3_500, :), 'tail', 'right');
                    p_lab1 = signrank(sust_data(idx_P1, :), sust_data(idx_P3_500, :), 'tail', 'left');
                    p_ab1  = min(p_rab1, p_lab1);
                    
                    if ~isnan(p_ab1) && p_ab1 < 0.05
                        x1 = get_x_pos(b_idx, idx_P1); x2 = get_x_pos(b_idx, idx_P3_500);
                        y = max(Y_mean(b_idx, [idx_P1, idx_P3_500]) + Y_se(b_idx, [idx_P1, idx_P3_500])) + y_shift;
                        plot([x1, x1, x2, x2], [y, y+1, y+1, y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                        star = '*'; if p_ab1 < 0.01, star = '**'; end; if p_ab1 < 0.001, star = '***'; end
                        text(mean([x1, x2]), y+2, star, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                    end
                end
            end
        end
    end
    
    % 3. Beta/Alpha Ratio (b=4): Unpred 2000 vs Cued & Unpred 2000 vs Unpred 500
    if isfield(aggregate_sust, 'BetaAlphaRatio') && length(aggregate_sust.BetaAlphaRatio.cohort) >= g && ~isempty(idx_P2_2000)
        ratio_data = aggregate_sust.BetaAlphaRatio.cohort(g).subj_data;
        if size(ratio_data, 2) > 1
            y_shift = 4;
            
            if ~isempty(idx_P1)
                p_r3 = signrank(ratio_data(idx_P2_2000, :), ratio_data(idx_P1, :), 'tail', 'right');
                p_l3 = signrank(ratio_data(idx_P2_2000, :), ratio_data(idx_P1, :), 'tail', 'left');
                p3   = min(p_r3, p_l3);
                
                if ~isnan(p3) && p3 < 0.05
                    x1 = get_x_pos(4, idx_P1); x2 = get_x_pos(4, idx_P2_2000);
                    y = max(Y_mean(4, [idx_P1, idx_P2_2000]) + Y_se(4, [idx_P1, idx_P2_2000])) + y_shift;
                    plot([x1, x1, x2, x2], [y, y+1.5, y+1.5, y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    star = '*'; if p3 < 0.01, star = '**'; end; if p3 < 0.001, star = '***'; end
                    text(mean([x1, x2]), y+3, star, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                    y_shift = y_shift + 8;
                end
            end
            
            if ~isempty(idx_P2_500)
                p_r4 = signrank(ratio_data(idx_P2_2000, :), ratio_data(idx_P2_500, :), 'tail', 'right');
                p_l4 = signrank(ratio_data(idx_P2_2000, :), ratio_data(idx_P2_500, :), 'tail', 'left');
                p4   = min(p_r4, p_l4);
                
                if ~isnan(p4) && p4 < 0.05
                    x1 = get_x_pos(4, idx_P2_500); x2 = get_x_pos(4, idx_P2_2000);
                    y = max(Y_mean(4, [idx_P2_500, idx_P2_2000]) + Y_se(4, [idx_P2_500, idx_P2_2000])) + y_shift;
                    plot([x1, x1, x2, x2], [y, y+1.5, y+1.5, y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    star = '*'; if p4 < 0.01, star = '**'; end; if p4 < 0.001, star = '***'; end
                    text(mean([x1, x2]), y+3, star, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                end
            end
        end
    end
    
    save_file_5d = fullfile(output_path, sprintf('Figure_5D_MultiBand_Sustained_%s', fig_tags{g}));
    try pause(0.5); saveas(fig5D, [save_file_5d, '.svg']); saveas(fig5D, [save_file_5d, '.png']); catch, end
    close(fig5D);
end
disp('All analyses completed successfully (All Retained Subjects, excluding S3 & S14).');

%% =========================================================================
%% GENERATE PNAS TABLE 1: SUMMARY OF GAIN, GEOMETRY, TIMING & SUSTAINED SHARE
%% =========================================================================
fprintf('\n=========================================================================================================\n');
fprintf(' PNAS TABLE 1: FREQUENCY-SPECIFIC DIVISION OF LABOR ACROSS GAIN, GEOMETRY, TIMING & SUSTAINED SHARE\n');
fprintf('=========================================================================================================\n');

t1_band = {}; t1_cond = {}; t1_gain = {}; t1_geom = {}; t1_lat = {}; t1_sust = {};
for b = 1:length(plot_bands)
    b_name = plot_bands{b};
    G_mat  = aggregate_G.(b_name).cohort(1).subj_data;    % [5 x n_ret]
    LR_mat = aggregate_log_r.(b_name).cohort(1).subj_data;% [5 x m x n_ret]
    LT_mat = aggregate_lat.(b_name).cohort(1).subj_data;  % [5 x n_ret]
    SU_mat = aggregate_sust.(b_name).cohort(1).subj_data; % [5 x n_ret]
    
    for c_i = 1:length(cfg.all_conds)
        c_str = cfg.all_conds{c_i};
        
        % 1. Gain (G +/- SEM, one-tailed signrank p vs P1 matching Figure 4B)
        g_vals = G_mat(c_i, :);
        p_g    = 1.000; 
        if c_i > 1
            p_g_r = signrank(g_vals, G_mat(1, :), 'tail', 'right');
            p_g_l = signrank(g_vals, G_mat(1, :), 'tail', 'left');
            p_g   = min(p_g_r, p_g_l);
        end
        g_str = sprintf('%.2f +/- %.2f (p=%.3f)', mean(g_vals,'omitnan'), std(g_vals,0,'omitnan')/sqrt(length(g_vals)), p_g);
        
        % 2. Geometry range across Xi_1..Xi_m
        lr_means = mean(squeeze(LR_mat(c_i, 1:cfg.m, :)), 2, 'omitnan');
        geom_str = sprintf('[%+.2f, %+.2f]', min(lr_means), max(lr_means));
        
        % 3. Latency shift (only for stimulus-delivered P1, P2_500, P3_500, one-tailed)
        if ismember(c_str, {'P1', 'P2_500', 'P3_500'})
            lt_vals = LT_mat(c_i, :);
            p_lt    = 1.000; 
            if c_i > 1
                p_lt_r = signrank(lt_vals, 0, 'tail', 'right');
                p_lt_l = signrank(lt_vals, 0, 'tail', 'left');
                p_lt   = min(p_lt_r, p_lt_l);
            end
            lat_str = sprintf('%+.1f +/- %.1f ms (p=%.3f)', mean(lt_vals,'omitnan'), std(lt_vals,0,'omitnan')/sqrt(length(lt_vals)), p_lt);
        else
            lat_str = 'N/A (Withheld)';
        end
        
        % 4. Sustained Share (%, one-tailed signrank p vs P1)
        su_vals = SU_mat(c_i, :);
        p_su    = 1.000; 
        if c_i > 1
            p_su_r = signrank(su_vals, SU_mat(1, :), 'tail', 'right');
            p_su_l = signrank(su_vals, SU_mat(1, :), 'tail', 'left');
            p_su   = min(p_su_r, p_su_l);
        end
        sust_str = sprintf('%.2f +/- %.2f%% (p=%.3f)', mean(su_vals,'omitnan'), std(su_vals,0,'omitnan')/sqrt(length(su_vals)), p_su);
        
        t1_band{end+1,1} = b_name; t1_cond{end+1,1} = c_str;
        t1_gain{end+1,1} = g_str;  t1_geom{end+1,1} = geom_str;
        t1_lat{end+1,1}  = lat_str; t1_sust{end+1,1} = sust_str;
        
        fprintf('%-14s | %-10s | Gain: %-24s | Geom: %-14s | Lat: %-26s | Sust: %s\n', ...
            b_name, c_str, g_str, geom_str, lat_str, sust_str);
    end
end
Table1_PNAS = table(t1_band, t1_cond, t1_gain, t1_geom, t1_lat, t1_sust, ...
    'VariableNames', {'Band', 'Condition', 'Global_Gain_G', 'Geometry_LogR_Range', 'Delta_Latency_ms', 'Sustained_Share_Pct'});
writetable(Table1_PNAS, fullfile(output_path, 'Table1_PNAS_Division_of_Labor.csv'));

%% =========================================================================
%% BAYES FACTOR (BF01) & TOST EQUIVALENCE CHECK FOR BETA GAIN & GEOMETRY
%% =========================================================================
fprintf('\n--- BAYES FACTOR (BF_01, Evidence for Null Stability) IN BETA BAND ---\n');
G_beta_subjs = aggregate_G.Beta.cohort(1).subj_data; % [5 conds x 12 subjs]
n_s = size(G_beta_subjs, 2);
for c_i = 2:length(cfg.all_conds)
    diff_G = G_beta_subjs(c_i, :) - G_beta_subjs(1, :);
    [~, p_t, ~, stats_t] = ttest(diff_G, 0);
    p_sr = signrank(G_beta_subjs(c_i, :), G_beta_subjs(1, :));
    %BIC / JZS approximation for BF01 (Rouder et al., 2009; Wagenmakers, 2007)
    bf01_G = sqrt(n_s) * (1 + (stats_t.tstat^2) / (n_s - 1))^(-n_s / 2);
    fprintf('Beta Gain %-12s vs P1: Mean G = %.3f | signrank p = %.3f | t = %+.2f | BF_01 = %.2f\n', ...
        cfg.all_conds{c_i}, mean(G_beta_subjs(c_i, :)), p_sr, stats_t.tstat, bf01_G);
end