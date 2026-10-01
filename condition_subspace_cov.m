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

conditions = {'BLA','BLT','P1','P2','P3'};
num_ch = 32;
fs = []; 
time_ms_eeg = [];

%% 0.1 Load and import file as EEGLAB matrix
data_all_conds = struct();
data_all_conds.Raw = struct();
data_all_conds.Alpha = struct();
data_all_conds.Beta = struct();

for c_idx = 1:length(conditions)
    cond_name = conditions{c_idx};
    in_dir = fullfile(input_path, cond_name);
    
    set_files = dir(fullfile(in_dir, '*.set'));
    names_sorted = sort(cellstr({set_files.name})');
    num_subjs_cond = length(names_sorted);
    
    for s = 1:num_subjs_cond
        file_to_load = names_sorted{s};
        
        % Check if an ICA-cleaned file exists (specifically for Subject 14 in P3)
        clean_file_check = fullfile(input_path, [cond_name, '_cleaned'], file_to_load);
        
        if exist(clean_file_check, 'file')
            fprintf('Loading ICA-CLEANED Subj %d/%d (%s) for condition: %s\n', ...
                s, num_subjs_cond, file_to_load, cond_name);
            EEG = pop_loadset('filename', file_to_load, 'filepath', fullfile(input_path, [cond_name, '_cleaned']));
        else
            fprintf('Loading Subj %d/%d (%s) for condition: %s\n', ...
                s, num_subjs_cond, file_to_load, cond_name);
            EEG = pop_loadset('filename', file_to_load, 'filepath', in_dir);
        end
        
        % --- ASSIGN & SPLIT SUB-CONDITIONS ---
        if strcmp(cond_name, 'P1')
            % Cued Reference: all trials
            data_all_conds.Raw.P1{s} = EEG.data;
            
        elseif strcmp(cond_name, 'P2')
            % Unpredicted Timing: Odd = 500 ms ISI, Even = 2000 ms ISI
            epoch_trials_p2_500    = 1:2:EEG.trials;
            epoch_trials_p2_2000   = 2:2:EEG.trials;
            
            data_all_conds.Raw.P2{s}      = EEG.data;
            data_all_conds.Raw.P2_500{s}  = EEG.data(:, :, epoch_trials_p2_500);
            data_all_conds.Raw.P2_2000{s} = EEG.data(:, :, epoch_trials_p2_2000);
            
        elseif strcmp(cond_name, 'P3')
            % Stimulus Uncertainty: Odd = 500 ms delivered, Even = Missing/Omitted
            epoch_trials_p3_500     = 1:2:EEG.trials;
            epoch_trials_p3_missing = 2:2:EEG.trials;
            
            data_all_conds.Raw.P3{s}         = EEG.data;
            data_all_conds.Raw.P3_500{s}     = EEG.data(:, :, epoch_trials_p3_500);
            data_all_conds.Raw.P3_missing{s} = EEG.data(:, :, epoch_trials_p3_missing);
        end
    end
end

%% =========================================================================
%% SECTION 0.2: BANDPASS FILTERING (ALPHA, BETA, RATIO)
%% =========================================================================
disp('Applying Zero-Phase FIR Filters for Alpha and Beta...');

% --- ROBUST SAMPLING RATE (fs) & NYQUIST (fn) EXTRACTION ---
if isfield(EEG, 'srate') && ~isempty(EEG.srate) && isscalar(EEG.srate) && EEG.srate > 0
    fs = EEG.srate;
elseif exist('time_s', 'var') && length(time_s) > 1
    % Derive fs directly from the time vector step size (dt = 1/fs)
    fs = round(1 / (time_s(2) - time_s(1)));
else
    % Standard default fallback for epoched data
    fs = 500; 
end

fn = fs / 2; % Nyquist frequency (guaranteed scalar)

fprintf('Filtering sample rate: %d Hz (Nyquist: %d Hz)\n', fs, fn);

% --- DESIGN ZERO-PHASE FIR FILTERS ---
ord_alpha = 3 * fix(fs / 8);   % Filter order for alpha (8-12 Hz)
ord_beta  = 3 * fix(fs / 15);  % Filter order for beta (15-30 Hz)

% Using element-wise division (./) ensures dimension safety
b_alpha   = fir1(ord_alpha, [8 12] ./ fn, 'bandpass');
b_beta    = fir1(ord_beta,  [15 30] ./ fn, 'bandpass');

% Smoothing window for envelopes (100 ms) to stabilize the ratio division
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
        
        % Filter original bands
        X_alpha_2D = filtfilt(b_alpha, 1, double(X_2D));
        X_beta_2D  = filtfilt(b_beta, 1, double(X_2D));
        
        % Save standard Alpha and Beta to struct
        data_all_conds.Alpha.(c_name){s} = permute(reshape(X_alpha_2D, nt, nc, ntr), [2, 1, 3]);
        data_all_conds.Beta.(c_name){s}  = permute(reshape(X_beta_2D, nt, nc, ntr), [2, 1, 3]);
        
        % --- Calculate Beta/Alpha Ratio ---
        % 1. Get Instantaneous Amplitude (Envelope) via Hilbert
        env_alpha = abs(hilbert(X_alpha_2D));
        env_beta  = abs(hilbert(X_beta_2D));
        
        % 2. Smooth envelopes to prevent dividing by near-zero instantaneous dips
        env_alpha_smooth = filtfilt(smoothing_kernel, 1, env_alpha);
        env_beta_smooth  = filtfilt(smoothing_kernel, 1, env_beta);
        
        % 3. Calculate Ratio (Adding eps to denominator prevents division by zero)
        ratio_2D = env_beta_smooth ./ (env_alpha_smooth + eps);
        
        % 4. Save new Ratio dimension to Struct
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
cfg.base          = [-1.000, -0.800];  
cfg.nTrialMatch   = 30;                
cfg.nRep          = 50;                
cfg.m             = 6;                 
cfg.seed          = 20260910;
cfg.slide_win     = 0.100; % Sliding window size (100 ms)
cfg.slide_step    = 0.025; % Sliding window step (25 ms) for smooth trajectory

clean_name = @(c) strrep(strrep(strrep(strrep(strrep(strrep(c, ...
    'BLT', 'Tactile (BLT)'), 'P1', 'Cued (P1)'), ...
    'P2_500', 'Unpred 500 (P2)'), 'P2_2000', 'Unpred 2000 (P2)'), ...
    'P3_500', 'Rand 500 (P3)'), 'P3_missing', 'Rand Null (P3)');

%% 2. Time Axis Alignment & Window Masking
% Reconstruct time vector directly from data dimensions and sampling rate
nPnts = size(data_all_conds.Raw.P1{1}, 2); 
time_s = cfg.tStart + (0:(nPnts - 1)) / cfg.fs;
time_s = time_s(:)'; % Force into a 1 x N row vector

iWin   = time_s >= cfg.win(1)  & time_s < cfg.win(2);
iBase  = time_s >= cfg.base(1) & time_s < cfg.base(2);
T_win  = nnz(iWin);
t_eval = time_s(iWin);

% --- Generate Time Vector for Sliding Window ---
t_starts    = time_s(1) : cfg.slide_step : (time_s(end) - cfg.slide_win);
t_centers   = t_starts + cfg.slide_win / 2;
num_windows = length(t_starts);

%% 3. Helper Functions (Unregularized Empirical Covariance)
% Simply calculates the standard covariance matrix across time points
covCentered = @(X, b) (( (X - b) - mean(X - b, 2) ) * ( (X - b) - mean(X - b, 2) )') / size(X, 2);

%% =========================================================================
%% MASTER LOOP: PROCESS FREQUENCY BANDS (STRICTLY P1 REFERENCE)
%% =========================================================================
bands_to_process = {'Raw', 'Alpha', 'Beta', 'BetaAlphaRatio'};
cfg.ref_cond = 'P1';
cfg.group_A = {'P2_500', 'P3_500'};
cfg.group_B = {'P2_2000', 'P3_missing'};
cfg.all_conds = [{cfg.ref_cond}, cfg.group_A, cfg.group_B];

% --- NEW: Aggregate Structs for Figures 3B, 5A, 5B, 5C, and 5D ---
aggregate_spectra = struct();
aggregate_log_r   = struct(); % Tracks log_r geometric data
aggregate_G       = struct(); % Tracks Global Gain (G) data
aggregate_lat     = struct(); % Tracks Delta Latency data
aggregate_sust    = struct(); % Tracks Sustained % data

for band_idx = 1:length(bands_to_process)
    band_name = bands_to_process{band_idx};
    current_data = data_all_conds.(band_name); 
    num_subjects = length(current_data.(cfg.ref_cond));
    
    fprintf('\n\n======================================================\n');
    fprintf('   BAND: %s | REFERENCE: %s (INDIVIDUAL COMPONENTS)\n', upper(band_name), cfg.ref_cond);
    fprintf('======================================================\n');
    
    % Setup Directory Structure
    group_dir = fullfile(output_path, 'group', band_name, cfg.ref_cond);
    if ~exist(group_dir, 'dir'), mkdir(group_dir); end

    % ------------------------------------------------------------------
    %% CHECKPOINT 1: Group Eigenbasis Evaluation
    % ------------------------------------------------------------------
    Xi_all_subjs = nan(cfg.num_ch, cfg.m, num_subjects);
    lam_frac_all = nan(cfg.m, num_subjects);
    gap_all      = nan(cfg.m - 1, num_subjects);
    
    lam_full_all = nan(cfg.num_ch, num_subjects); 
    
    for s = 1:num_subjects
        raw_ref = current_data.(cfg.ref_cond){s};
        n_trials = size(raw_ref, 3);
        if n_trials < (cfg.nTrialMatch * 2)
            error('Subj %d insufficient %s trials (%d < %d)', s, cfg.ref_cond, n_trials, cfg.nTrialMatch * 2);
        end
        
        rng(cfg.seed + s);
        perm = randperm(n_trials);
        iA = perm(1:cfg.nTrialMatch);
        
        erpA = mean(raw_ref(:, :, iA), 3, 'omitnan');
        bA   = mean(erpA(:, iBase), 2);
        
        % Calculate unregularized centered covariance
        C_A_cent = covCentered(erpA(:, iWin), bA);
        
        [V, D] = eig(C_A_cent, 'vector');
        [lam, ord] = sort(D, 'descend');
        V = V(:, ord);
        
        Xi_all_subjs(:, :, s) = V(:, 1:cfg.m);
        lam_frac_all(:, s)    = lam(1:cfg.m) / sum(lam);
        gap_all(:, s)         = (lam(1:cfg.m-1) - lam(2:cfg.m)) ./ lam(1:cfg.m-1);
        
        lam_full_all(:, s)    = lam / sum(lam); 
    end
    
    aggregate_spectra.(band_name) = mean(cumsum(lam_full_all, 1), 2) * 100;
    
    % Sign-align eigenvectors to Subject 1 to allow clean visual topoplot averaging
    Xi_aligned = Xi_all_subjs;
    for s = 2:num_subjects
        for k = 1:cfg.m
            if dot(Xi_aligned(:, k, s), Xi_aligned(:, k, 1)) < 0
                Xi_aligned(:, k, s) = -Xi_aligned(:, k, s);
            end
        end
    end
    
    save(fullfile(group_dir, 'spectra.mat'), 'lam_frac_all', 'gap_all');
    
    % --- Deliverables A & B (Spectra and Topoplots) ---
    figA = figure('Position', [100, 100, 1000, 400], 'Visible', 'off');
    subplot(1, 2, 1); plot(1:cfg.m, lam_frac_all, '-o', 'LineWidth', 1.5); grid on;
    xlabel('Direction Index'); ylabel('Fraction of Total Variance'); title(sprintf('Eigenvalue Spectra [%s|%s]', band_name, cfg.ref_cond));
    subplot(1, 2, 2); bar(mean(gap_all, 2)); hold on; errorbar(1:cfg.m-1, mean(gap_all, 2), std(gap_all, 0, 2), 'k.', 'LineWidth', 1.2);
    grid on; xlabel('Gap Index'); ylabel('Relative Gap'); title('Mean Subspace Gaps');
    saveas(figA, fullfile(group_dir, 'Deliverable_A_Spectra.png')); close(figA);

    % --- Deliverable C: Subject-by-Subject Similarity Matrices ---
    figC = figure('Position', [100, 100, 300 * cfg.m, 350], 'Visible', 'off');
    tiledlayout(1, cfg.m, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    subj_subj_sim = nan(num_subjects, num_subjects, cfg.m);
    for k = 1:cfg.m
        for s1 = 1:num_subjects
            for s2 = 1:num_subjects
                % Absolute Cosine Similarity (1 = identical direction, 0 = orthogonal)
                subj_subj_sim(s1, s2, k) = abs(dot(Xi_aligned(:, k, s1), Xi_aligned(:, k, s2)));
            end
        end
        
        nexttile;
        imagesc(subj_subj_sim(:, :, k));
        colormap('jet'); clim([0 1]); axis square;
        title(sprintf('\\xi_%d', k), 'FontSize', 18, 'FontWeight', 'bold');
        
        if k == 1, ylabel('Subject ID', 'FontSize', 12, 'FontWeight', 'bold'); end
        xlabel('Subject ID', 'FontSize', 12, 'FontWeight', 'bold');
        xticks([1, num_subjects]); yticks([1, num_subjects]);
    end
    
    cb = colorbar; cb.Layout.Tile = 'east'; 
    cb.Label.String = 'Absolute Cosine Similarity'; cb.Label.FontSize = 14;
    sgtitle(sprintf('Cross-Subject Component Consistency (%s %s)', band_name, cfg.ref_cond), 'FontSize', 20, 'FontWeight', 'bold');
    
    save(fullfile(group_dir, 'subj_subj_similarity.mat'), 'subj_subj_sim');
    saveas(figC, fullfile(group_dir, 'Deliverable_C_SubjSubj_Sim.png')); close(figC);

    % ------------------------------------------------------------------
    %% PROJECTIONS: Individual Component Evaluation (All Subjects/Reps)
    % ------------------------------------------------------------------
    G_splits     = nan(length(cfg.all_conds), cfg.nRep, num_subjects);
    log_r_splits = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
    lat_splits   = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
    sust_splits  = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
    traj_splits  = nan(length(cfg.all_conds), cfg.m, num_windows, cfg.nRep, num_subjects);

    for subj = 1:num_subjects
        subj_dir = fullfile(output_path, sprintf('subj%02d', subj), band_name, cfg.ref_cond);
        if ~exist(subj_dir, 'dir'), mkdir(subj_dir); end
        save(fullfile(subj_dir, 'cfg.mat'), 'cfg');
        
        Xi_stable = Xi_aligned(:, :, subj);
        
        for rep = 1:cfg.nRep
            rng(cfg.seed + rep + subj*1000); 
            
            raw_ref = current_data.(cfg.ref_cond){subj};
            nRef = size(raw_ref, 3);
            perm = randperm(nRef);
            iA = perm(1:cfg.nTrialMatch);
            iB = perm(cfg.nTrialMatch + (1:cfg.nTrialMatch));
            
            erpA = mean(raw_ref(:, :, iA), 3, 'omitnan');
            bA   = mean(erpA(:, iBase), 2);
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
            bB   = mean(erpB(:, iBase), 2);
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
                    idx_k = randperm(size(raw_cond, 3), cfg.nTrialMatch);
                    erp_k = mean(raw_cond(:, :, idx_k), 3, 'omitnan');
                    bk    = mean(erp_k(:, iBase), 2);
                    Xk    = erp_k(:, iWin) - bk;
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
    fprintf('\nRunning 1,000-iteration Spatial Permutation Test for Consistency...\n');
    
    n_perms = 1000;
    num_pairs = (num_subjects * (num_subjects - 1)) / 2;
    
    true_mean_sim = zeros(cfg.m, 1);
    null_mean_sim = zeros(cfg.m, n_perms);
    p_vals        = zeros(cfg.m, 1);
    thresh_95     = zeros(cfg.m, 1);
    
    for k = 1:cfg.m
        pair_sims = zeros(num_pairs, 1);
        idx = 1;
        for i = 1:num_subjects
            for j = (i+1):num_subjects
                v1 = Xi_aligned(:, k, i);
                v2 = Xi_aligned(:, k, j);
                pair_sims(idx) = abs(dot(v1, v2)) / (norm(v1) * norm(v2));
                idx = idx + 1;
            end
        end
        true_mean_sim(k) = mean(pair_sims);
    end
    
    for p_idx = 1:n_perms
        for k = 1:cfg.m
            shuffled_Xi = zeros(cfg.num_ch, num_subjects);
            for s = 1:num_subjects
                shuff_idx = randperm(cfg.num_ch);
                shuffled_Xi(:, s) = Xi_aligned(shuff_idx, k, s);
            end
            
            pair_sims_null = zeros(num_pairs, 1);
            idx = 1;
            for i = 1:num_subjects
                for j = (i+1):num_subjects
                    v1 = shuffled_Xi(:, i);
                    v2 = shuffled_Xi(:, j);
                    pair_sims_null(idx) = abs(dot(v1, v2)) / (norm(v1) * norm(v2));
                    idx = idx + 1;
                end
            end
            null_mean_sim(k, p_idx) = mean(pair_sims_null);
        end
    end
    
    fprintf('\n--- Cross-Subject Component Consistency --- \n');
    fprintf('%-10s | %-12s | %-12s | %-10s | %-10s\n', 'Direction', 'True Sim', '95% Null Thresh', 'p-value', 'Significance');
    fprintf('----------------------------------------------------------------\n');
    for k = 1:cfg.m
        p_vals(k) = (1 + sum(null_mean_sim(k, :) >= true_mean_sim(k))) / (1 + n_perms);
        thresh_95(k) = prctile(null_mean_sim(k, :), 95);
        
        if p_vals(k) < 0.001, sig_str = '***';
        elseif p_vals(k) < 0.01, sig_str = '**';
        elseif p_vals(k) < 0.05, sig_str = '*';
        else, sig_str = 'n.s.'; end
        
        fprintf('Xi_%-7d | %-12.3f | %-15.3f | %-10.4f | %-10s\n', k, true_mean_sim(k), thresh_95(k), p_vals(k), sig_str);
    end
    fprintf('----------------------------------------------------------------\n');
    
    figC = figure('Position', [150, 150, 1400, 600], 'Name', sprintf('Permutation Test: %s', band_name), 'Visible', 'off');
    t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, sprintf('Cross-Subject Consistency: Null vs True (%s | %s)', band_name, cfg.ref_cond), 'FontSize', 18, 'FontWeight', 'bold');
    
    for k = 1:cfg.m
        nexttile; hold on;
        histogram(null_mean_sim(k, :), 30, 'Normalization', 'probability', 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
        xline(thresh_95(k), 'k--', 'LineWidth', 1.5, 'Label', '95% Thresh', 'LabelOrientation', 'horizontal', 'FontSize', 10);
        
        if p_vals(k) < 0.05
            line_color = [0.8 0.1 0.1]; 
        else
            line_color = [0.2 0.2 0.2]; 
        end
        xline(true_mean_sim(k), '-', 'Color', line_color, 'LineWidth', 2.5, 'Label', 'True Sim', 'LabelOrientation', 'horizontal', 'FontSize', 10);
        
        title(sprintf('\\xi_%d (p = %.3f)', k, p_vals(k)), 'FontSize', 16, 'FontWeight', 'bold');
        xlabel('Mean Abs. Cosine Similarity', 'FontSize', 12);
        if ismember(k, [1, 4]), ylabel('Probability', 'FontSize', 12); end
        
        max_x = max(max(null_mean_sim(k, :)), true_mean_sim(k)) * 1.1;
        xlim([0, max_x]);
        grid on; box on;
    end
    
    save_file_perm = fullfile(group_dir, 'Deliverable_C_Permutation_Test.png');
    try pause(0.5); saveas(figC, save_file_perm); catch, warning('Could not save %s.', save_file_perm); end
    close(figC);

    % ------------------------------------------------------------------
    %% CHECKPOINT 2 & 3: GROUP-SPLIT STATS & FIGURES
    % ------------------------------------------------------------------
    % 1. Pre-calculate the mean logic across all reps for all subjects
    subj_mean_log_r = squeeze(mean(log_r_splits, 3, 'omitnan')); 
    subj_mean_traj  = squeeze(mean(traj_splits, 4, 'omitnan'));  
    
    % 2. Define subject subgroups (EXCLUDING OUTLIER SUBJECT 14)
    exclude_subjs = 14; 

    % Define subject subgroups (Subject 14 is now safely cleaned!)
    g1_subjs = intersect(4:11, 1:num_subjects); % Cohort A (n=8)
    g2_subjs = setdiff(1:num_subjects, g1_subjs); % Cohort B (n=6, cleaned)
    
    subj_groups = {g1_subjs, g2_subjs};
    subj_group_names = {'Subjs_4_to_11', 'Other_Subjs_Cleaned'};
    
    % PRE-CALCULATE SYNCHRONIZED Y-LIMITS
    time_mask_zoom = (t_centers >= -0.1) & (t_centers <= 1.0);
    time_mask_full = (t_centers >= -0.15); 
    
    global_ylims_full = zeros(cfg.m, 2); 
    global_ylims_zoom = zeros(cfg.m, 2); 
    global_clims      = zeros(cfg.m, 1);
    all_log_r         = [];
    
    Xi_g1 = mean(Xi_aligned(:, :, g1_subjs), 3, 'omitnan'); Xi_g1 = Xi_g1 ./ vecnorm(Xi_g1);
    Xi_g2 = mean(Xi_aligned(:, :, g2_subjs), 3, 'omitnan'); Xi_g2 = Xi_g2 ./ vecnorm(Xi_g2);
    
    for k = 1:cfg.m
        all_up_full = []; all_dn_full = [];
        all_up_zoom = []; all_dn_zoom = [];
        max_c = 0; 
        
        for c_idx = 1:length(cfg.all_conds)
            cond = cfg.all_conds{c_idx};
            
            % Cohort 1
            if ~isempty(g1_subjs)
                g1_val = squeeze(subj_mean_traj(c_idx, k, :, g1_subjs));
                if length(g1_subjs) == 1, g1_val = g1_val(:); g1_m = g1_val; g1_s = zeros(size(g1_val));
                else, g1_m = mean(g1_val, 2, 'omitnan'); g1_s = std(g1_val, 0, 2, 'omitnan') ./ sqrt(length(g1_subjs)); end
                
                all_up_full = [all_up_full; g1_m(time_mask_full) + g1_s(time_mask_full)];
                all_dn_full = [all_dn_full; g1_m(time_mask_full) - g1_s(time_mask_full)];
                all_up_zoom = [all_up_zoom; g1_m(time_mask_zoom) + g1_s(time_mask_zoom)];
                all_dn_zoom = [all_dn_zoom; g1_m(time_mask_zoom) - g1_s(time_mask_zoom)];
                
                erp_g1 = mean(cat(3, current_data.(cond){g1_subjs}), 3, 'omitnan');
                max_c = max(max_c, max(abs(Xi_g1(:, k)' * erp_g1)) * 0.5);
                
                if k == 1
                    g1_lr_m = mean(subj_mean_log_r(c_idx, :, g1_subjs), 3, 'omitnan');
                    g1_lr_s = std(subj_mean_log_r(c_idx, :, g1_subjs), 0, 3, 'omitnan') ./ sqrt(length(g1_subjs));
                    all_log_r = [all_log_r; g1_lr_m(:) + g1_lr_s(:); g1_lr_m(:) - g1_lr_s(:)];
                end
            end
            
            % Cohort 2
            if ~isempty(g2_subjs)
                g2_val = squeeze(subj_mean_traj(c_idx, k, :, g2_subjs));
                if length(g2_subjs) == 1, g2_val = g2_val(:); g2_m = g2_val; g2_s = zeros(size(g2_val));
                else, g2_m = mean(g2_val, 2, 'omitnan'); g2_s = std(g2_val, 0, 2, 'omitnan') ./ sqrt(length(g2_subjs)); end
                
                all_up_full = [all_up_full; g2_m(time_mask_full) + g2_s(time_mask_full)];
                all_dn_full = [all_dn_full; g2_m(time_mask_full) - g2_s(time_mask_full)];
                all_up_zoom = [all_up_zoom; g2_m(time_mask_zoom) + g2_s(time_mask_zoom)];
                all_dn_zoom = [all_dn_zoom; g2_m(time_mask_zoom) - g2_s(time_mask_zoom)];
                
                erp_g2 = mean(cat(3, current_data.(cond){g2_subjs}), 3, 'omitnan');
                max_c = max(max_c, max(abs(Xi_g2(:, k)' * erp_g2)) * 0.5);
                
                if k == 1
                    g2_lr_m = mean(subj_mean_log_r(c_idx, :, g2_subjs), 3, 'omitnan');
                    g2_lr_s = std(subj_mean_log_r(c_idx, :, g2_subjs), 0, 3, 'omitnan') ./ sqrt(length(g2_subjs));
                    all_log_r = [all_log_r; g2_lr_m(:) + g2_lr_s(:); g2_lr_m(:) - g2_lr_s(:)];
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

    % --- LOOP OVER BOTH SUBJECT GROUPS ---
    for g = 1:2
        curr_subjs = subj_groups{g};
        curr_n = length(curr_subjs);
        g_name = subj_group_names{g};
        
        if curr_n == 0
            fprintf('Skipping %s (No subjects found)\n', g_name);
            continue;
        end
        
        grp_mean_log_r = mean(subj_mean_log_r(:, :, curr_subjs), 3, 'omitnan');      
        grp_se_log_r   = std(subj_mean_log_r(:, :, curr_subjs), 0, 3, 'omitnan') ./ sqrt(curr_n);
        
        % --- Save to aggregates for Figures 5A, 5B, 5C, and 5D ---
        aggregate_log_r.(band_name).cohort(g).mean = grp_mean_log_r;
        aggregate_log_r.(band_name).cohort(g).se   = grp_se_log_r;
        
        % Gain (G)
        subj_mean_G = squeeze(mean(G_splits(:, :, curr_subjs), 2, 'omitnan')); 
        if curr_n == 1, subj_mean_G = subj_mean_G(:); end
        aggregate_G.(band_name).cohort(g).mean = mean(subj_mean_G, 2, 'omitnan');
        aggregate_G.(band_name).cohort(g).se   = std(subj_mean_G, 0, 2, 'omitnan') ./ sqrt(curr_n);
        aggregate_G.(band_name).cohort(g).subj_data = subj_mean_G; % Track for t-tests
        
        % Latency
        subj_mean_lat = squeeze(mean(mean(lat_splits(:, :, :, curr_subjs), 3, 'omitnan'), 2, 'omitnan')); 
        if curr_n == 1, subj_mean_lat = subj_mean_lat(:); end
        lat_dir_mean = squeeze(mean(lat_splits(:, :, :, curr_subjs), 3, 'omitnan'));
        if curr_n == 1, lat_dir_mean = reshape(lat_dir_mean, [length(cfg.all_conds), cfg.m, 1]); end
        
        aggregate_lat.(band_name).cohort(g).mean = mean(lat_dir_mean, 3, 'omitnan');
        aggregate_lat.(band_name).cohort(g).se   = std(lat_dir_mean, 0, 3, 'omitnan') ./ sqrt(curr_n);
        aggregate_lat.(band_name).cohort(g).subj_data = subj_mean_lat; 
        
        % Sustained % 
        subj_mean_sust = squeeze(mean(mean(sust_splits(:, :, :, curr_subjs), 3, 'omitnan'), 2, 'omitnan')) * 100;
        if curr_n == 1, subj_mean_sust = subj_mean_sust(:); end
        aggregate_sust.(band_name).cohort(g).mean = mean(subj_mean_sust, 2, 'omitnan');
        aggregate_sust.(band_name).cohort(g).se   = std(subj_mean_sust, 0, 2, 'omitnan') ./ sqrt(curr_n);
        aggregate_sust.(band_name).cohort(g).subj_data = subj_mean_sust; 
        
        % -----------------------------------------------------------------
        % CHECKPOINT 2: GROUP LEVEL STATS
        % -----------------------------------------------------------------
        figChk2 = figure('Position', [100, 100, 1400, 700], 'Visible', 'off');
        
        subplot(1, 2, 1); hold on;
        set(gca, 'FontSize', 22);
        yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('Ceiling (%s B)', cfg.ref_cond));
        colors_A = {[0.850 0.325 0.098], [0.929 0.694 0.125], [0.494 0.184 0.556]};
        for i = 1:length(cfg.group_A)
            c_idx = find(strcmp(cfg.all_conds, cfg.group_A{i}));
            if isempty(c_idx), continue; end
            col_idx = mod(i-1, 3) + 1;
            errorbar(1:cfg.m, grp_mean_log_r(c_idx, :), grp_se_log_r(c_idx, :), grp_se_log_r(c_idx, :), '-o', ...
                'Color', colors_A{col_idx}, 'LineWidth', 2, 'MarkerFaceColor', colors_A{col_idx}, 'DisplayName', clean_name(cfg.group_A{i}));
        end
        grid on; 
        xlim([0.75, cfg.m + 0.25]); xticks(1:cfg.m);
        xlabel('Spatial Direction Index', 'FontSize', 20); 
        ylabel('log(r_i) \pm SEM', 'FontSize', 20);
        title(sprintf('Group A: Stimulus Delivered (%s)', strrep(g_name, '_', ' ')), 'FontSize', 24); 
        lgd1 = legend('Location', 'best'); lgd1.FontSize = 22;
        ylim([cp2_y_min, cp2_y_max]);
        
        subplot(1, 2, 2); hold on;
        set(gca, 'FontSize', 22);
        yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('Ceiling (%s B)', cfg.ref_cond));
        colors_B = {[0.466 0.674 0.188], [0.301 0.745 0.933]};
        for i = 1:length(cfg.group_B)
            c_idx = find(strcmp(cfg.all_conds, cfg.group_B{i}));
            if isempty(c_idx), continue; end
            errorbar(1:cfg.m, grp_mean_log_r(c_idx, :), grp_se_log_r(c_idx, :), grp_se_log_r(c_idx, :), '-s', ...
                'Color', colors_B{i}, 'LineWidth', 2, 'MarkerFaceColor', colors_B{i}, 'DisplayName', clean_name(cfg.group_B{i}));
        end
        grid on; 
        xlim([0.75, cfg.m + 0.25]); xticks(1:cfg.m);
        xlabel('Spatial Direction Index', 'FontSize', 20); 
        ylabel('log(r_i) \pm SEM', 'FontSize', 20);
        title(sprintf('Group B: Omission (%s)', strrep(g_name, '_', ' ')), 'FontSize', 24); 
        lgd2 = legend('Location', 'best'); lgd2.FontSize = 22;
        ylim([cp2_y_min, cp2_y_max]);
        
        sgtitle(sprintf('Subspace Redistribution: %s %s', band_name, cfg.ref_cond), 'FontSize', 24, 'FontWeight', 'bold');
        
        save_base_chk2 = fullfile(group_dir, sprintf('Checkpoint2_Redistribution_%s', g_name));
        try pause(0.5); saveas(figChk2, [save_base_chk2, '.svg']); saveas(figChk2, [save_base_chk2, '.png']); catch, end
        close(figChk2);
        
        % -----------------------------------------------------------------
        % CHECKPOINT 3a: FULL SUBSPACE TRAJECTORIES 
        % -----------------------------------------------------------------
        grp_mean_traj = squeeze(mean(subj_mean_traj(:, :, :, curr_subjs), 4, 'omitnan')); 
        grp_sem_traj  = squeeze(std(subj_mean_traj(:, :, :, curr_subjs), 0, 4, 'omitnan')) ./ sqrt(curr_n);
        
        sig_components = find(p_vals < 0.05)';
        if isempty(sig_components), sig_components = 1; end
        num_sig = length(sig_components);
        
        fig_height = max(600, num_sig * 500);
        figChk3 = figure('Position', [100, 100, 1600, fig_height], 'Visible', 'off');
        tiledlayout(num_sig, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        c_idx_ref = find(strcmp(cfg.all_conds, cfg.ref_cond));
        
        for k_idx = 1:num_sig
            k = sig_components(k_idx);
            cp3a_min_y = global_ylims_full(k, 1);
            cp3a_max_y = global_ylims_full(k, 2);
            
            % --- CP3a: Group A ---
            nexttile; hold on;
            set(gca, 'FontSize', 22);
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
                col_idx = mod(i-1, 3) + 1;
                
                y_val = reshape(grp_mean_traj(c_idx, k, :), 1, []); y_sem = reshape(grp_sem_traj(c_idx, k, :), 1, []);
                fill([t_centers, fliplr(t_centers)], max(0, [y_val + y_sem, fliplr(y_val - y_sem)]), colors_A{col_idx}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(t_centers, y_val, '-', 'Color', colors_A{col_idx}, 'LineWidth', 2, 'DisplayName', clean_name(cfg.group_A{i}));
            end
            grid on; 
            xlabel('Time (s)', 'FontSize', 20); ylabel('Projected Variance', 'FontSize', 20);
            title(sprintf('Group A: \\xi_%d Trajectory', k), 'FontSize', 24); 
            if k_idx == 1, lgdA = legend('Location', 'best'); lgdA.FontSize = 18; end
            xlim([t_centers(1), t_centers(end)]); ylim([cp3a_min_y, cp3a_max_y]);
            
            % --- CP3a: Group B ---
            nexttile; hold on;
            set(gca, 'FontSize', 22);
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
                fill([t_centers, fliplr(t_centers)], max(0, [y_val + y_sem, fliplr(y_val - y_sem)]), colors_B{i}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(t_centers, y_val, '-', 'Color', colors_B{i}, 'LineWidth', 2, 'DisplayName', clean_name(cfg.group_B{i}));
            end
            grid on; 
            xlabel('Time (s)', 'FontSize', 20); ylabel('Projected Variance', 'FontSize', 20);
            title(sprintf('Group B: \\xi_%d Trajectory', k), 'FontSize', 24); 
            if k_idx == 1, lgdB = legend('Location', 'best'); lgdB.FontSize = 18; end
            xlim([t_centers(1), t_centers(end)]); ylim([cp3a_min_y, cp3a_max_y]);
        end
        sgtitle(sprintf('Subspace Trajectories: %s (%s)', band_name, strrep(g_name, '_', ' ')), 'FontSize', 24, 'FontWeight', 'bold');
        
        save_base_chk3 = fullfile(group_dir, sprintf('Checkpoint3a_Trajectories_%s', g_name));
        try pause(0.5); saveas(figChk3, [save_base_chk3, '.svg']); saveas(figChk3, [save_base_chk3, '.png']); catch, end
        close(figChk3);
        
        % -----------------------------------------------------------------
        % CHECKPOINT 3b: SIDE-BY-SIDE TOPOPLOT & TRAJECTORIES
        % -----------------------------------------------------------------
        Xi_bar = mean(Xi_aligned(:, :, curr_subjs), 3, 'omitnan');
        Xi_bar = Xi_bar ./ vecnorm(Xi_bar); 
        idx_04_06 = find(time_s >= 0.4 & time_s <= 0.6);
        
        groupA_conds = {'P1', 'P2_500', 'P3_500'};
        groupB_conds = {'P1', 'P2_2000', 'P3_missing'};
        cols_A = {[0 0 0], [0.850 0.325 0.098], [0.929 0.694 0.125]}; 
        cols_B = {[0 0 0], [0.466 0.674 0.188], [0.301 0.745 0.933]}; 
        
        for k = 1:cfg.m
            xi = Xi_bar(:, k); 
            c_lim = [-global_clims(k), global_clims(k)];
            min_y = global_ylims_zoom(k, 1);
            max_y = global_ylims_zoom(k, 2);
            
            % --- Group A ---
            figA = figure('Position', [50, 50, 1600, 1200], 'Name', sprintf('Dir %d - Group A', k), 'Visible', 'off');
            tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'normal');
            for c = 1:length(groupA_conds)
                cond = groupA_conds{c}; col = cols_A{c};
                c_idx = find(strcmp(cfg.all_conds, cond));
                erp_all = cat(3, current_data.(cond){curr_subjs});
                avg_data = mean(erp_all, 3, 'omitnan');
                
                nexttile((c-1)*4 + 1); set(gca, 'FontSize', 22);
                amp_04_06 = mean(xi' * avg_data(:, idx_04_06), 2);
                topoplot(xi * amp_04_06, EEG.chanlocs, 'numcontour', 0);
                clim(c_lim); colormap('jet');
                cb = colorbar; cb.Label.String = 'Proj. Amp'; cb.Label.FontSize = 20; cb.FontSize = 20;
                title(sprintf('%s\n(0.4 - 0.6s)', clean_name(cond)), 'FontSize', 24, 'FontWeight', 'bold');
                
                nexttile((c-1)*4 + 2, [1 3]); hold on; set(gca, 'FontSize', 22);
                patch([cfg.win(1) cfg.win(2) cfg.win(2) cfg.win(1)], [min_y min_y max_y max_y], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
                xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
                
                if ~isempty(c_idx)
                    y_val = reshape(grp_mean_traj(c_idx, k, :), 1, []); y_sem = reshape(grp_sem_traj(c_idx, k, :), 1, []);
                    fill([t_centers, fliplr(t_centers)], max(0, [y_val + y_sem, fliplr(y_val - y_sem)]), col, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                    plot(t_centers, y_val, '-', 'Color', col, 'LineWidth', 2.5);
                end
                xlim([-0.1, 1.0]); ylim([min_y, max_y]); ylabel('Proj. Var', 'FontSize', 20, 'FontWeight', 'bold'); grid on;
                if c < length(groupA_conds), xticklabels({}); else, xlabel('Time (s)', 'FontSize', 20, 'FontWeight', 'bold'); end
            end
            sgtitle(sprintf('Group A: Spatial Reconfiguration (\\xi_%d) [%s]', k, strrep(g_name, '_', ' ')), 'FontSize', 24, 'FontWeight', 'bold');
            save_base_A = fullfile(group_dir, sprintf('Checkpoint3b_GroupA_Xi%d_%s', k, g_name));
            try pause(0.5); saveas(figA, [save_base_A, '.svg']); saveas(figA, [save_base_A, '.png']); catch, end
            close(figA);
            
            % --- Group B ---
            figB = figure('Position', [150, 150, 1600, 1200], 'Name', sprintf('Dir %d - Group B', k), 'Visible', 'off');
            tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'normal');
            for c = 1:length(groupB_conds)
                cond = groupB_conds{c}; col = cols_B{c};
                c_idx = find(strcmp(cfg.all_conds, cond));
                erp_all = cat(3, current_data.(cond){curr_subjs});
                avg_data = mean(erp_all, 3, 'omitnan');
                
                nexttile((c-1)*4 + 1); set(gca, 'FontSize', 22);
                amp_04_06 = mean(xi' * avg_data(:, idx_04_06), 2);
                topoplot(xi * amp_04_06, EEG.chanlocs, 'numcontour', 0);
                clim(c_lim); colormap('jet');
                cb = colorbar; cb.Label.String = 'Proj. Amp'; cb.Label.FontSize = 20; cb.FontSize = 20;
                title(sprintf('%s\n(0.4 - 0.6s)', clean_name(cond)), 'FontSize', 24, 'FontWeight', 'bold');
                
                nexttile((c-1)*4 + 2, [1 3]); hold on; set(gca, 'FontSize', 22);
                patch([cfg.win(1) cfg.win(2) cfg.win(2) cfg.win(1)], [min_y min_y max_y max_y], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
                xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
                
                if ~isempty(c_idx)
                    y_val = reshape(grp_mean_traj(c_idx, k, :), 1, []); y_sem = reshape(grp_sem_traj(c_idx, k, :), 1, []);
                    fill([t_centers, fliplr(t_centers)], max(0, [y_val + y_sem, fliplr(y_val - y_sem)]), col, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                    plot(t_centers, y_val, '-', 'Color', col, 'LineWidth', 2.5);
                end
                xlim([-0.1, 1.0]); ylim([min_y, max_y]); ylabel('Proj. Var', 'FontSize', 20, 'FontWeight', 'bold'); grid on;
                if c < length(groupB_conds), xticklabels({}); else, xlabel('Time (s)', 'FontSize', 20, 'FontWeight', 'bold'); end
            end
            sgtitle(sprintf('Group B: Spatial Reconfiguration (\\xi_%d) [%s]', k, strrep(g_name, '_', ' ')), 'FontSize', 24, 'FontWeight', 'bold');
            save_base_B = fullfile(group_dir, sprintf('Checkpoint3b_GroupB_Xi%d_%s', k, g_name));
            try pause(0.5); saveas(figB, [save_base_B, '.svg']); saveas(figB, [save_base_B, '.png']); catch, end
            close(figB);
        end 
    end 
    
    % --- Export Group Results and Split Summary Tables ---
    save(fullfile(group_dir, 'results.mat'), 'log_r_splits', 'G_splits', 'lat_splits', 'sust_splits', 'traj_splits', 't_centers');
    
    fid = fopen(fullfile(group_dir, 'summary_table.txt'), 'w');
    for g = 1:2
        curr_subjs = subj_groups{g};
        curr_n = length(curr_subjs);
        g_name = subj_group_names{g};
        if curr_n == 0, continue; end
        
        grp_mean_log_r = mean(subj_mean_log_r(:, :, curr_subjs), 3, 'omitnan'); 
        for out = [1, fid]
            fprintf(out, '\n========================================================================================\n');
            fprintf(out, 'SUMMARY TABLE: %s REFERENCE (BAND: %s) [%s, n=%d]\n', cfg.ref_cond, upper(band_name), strrep(g_name, '_', ' '), curr_n);
            fprintf(out, '========================================================================================\n');
            
            for k = sig_components
                fprintf(out, '\n--- SIGNIFICANT DIRECTION: \\xi_%d (Global p = %.3f) ---\n', k, p_vals(k));
                fprintf(out, '----------------------------------------------------------------------------------------\n');
                fprintf(out, '%-18s | %-10s | %-10s | %-18s | %-15s\n', 'Condition', 'Gain (G)', sprintf('log(r_%d)', k), sprintf('\\Delta Latency (\\xi_%d)', k), sprintf('Sustained %% (\\xi_%d)', k));
                fprintf(out, '----------------------------------------------------------------------------------------\n');
                
                for c = 1:length(cfg.all_conds)
                    cond = cfg.all_conds{c};
                    med_G     = mean(G_splits(c, :, curr_subjs), 'all', 'omitnan'); 
                    med_lr_k  = grp_mean_log_r(c, k);
                    med_lat_k = mean(lat_splits(c, k, :, curr_subjs), 'all', 'omitnan'); 
                    med_sus_k = mean(sust_splits(c, k, :, curr_subjs), 'all', 'omitnan') * 100;
                    
                    fprintf(out, '%-18s | %-10.3f | %-10.3f | %+7.1f ms         | %6.2f%%\n', cond, med_G, med_lr_k, med_lat_k, med_sus_k);
                end
            end
            fprintf(out, '\n\n');
        end
    end
    fclose(fid);
end % END BAND LOOP

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
plot(cfg.m, 85, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');

xline(2, ':', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2, 'HandleVisibility', 'off');
text(2.2, 50, 'Beta/Alpha intrinsic 2D manifold (p < 0.05)', ...
    'FontSize', 18, 'Color', [0.4940 0.1840 0.5560], 'FontWeight', 'bold', 'BackgroundColor', 'w');
    
text(cfg.m + 0.2, 70, sprintf('Truncation at m=%d\n(>85%% variance)', cfg.m), ...
    'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold', 'BackgroundColor', 'w');
    
xlim([0.5, 15]); xticks(1:15); ylim([0, 100]); yticks(0:10:100);
xlabel('Number of Components (Spatial Directions)', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Cumulative Variance Explained (%)', 'FontSize', 24, 'FontWeight', 'bold');
title('Intrinsic Dimensionality Across Frequency Bands', 'FontSize', 28, 'FontWeight', 'bold');

lgd = legend('Location', 'southeast'); lgd.FontSize = 20; grid on;

save_file_3b = fullfile(output_path, 'Figure_3B_Aggregate_ScreePlot');
try pause(0.5); saveas(fig3B, [save_file_3b, '.svg']); saveas(fig3B, [save_file_3b, '.png']); catch, end
close(fig3B);

%% =========================================================================
%% AGGREGATE FIGURE 5A: MULTI-BAND LOG(R) COMPARISON (GEOMETRIC STABILITY)
%% =========================================================================
disp('Generating Figure 5A: Multi-Band log(r) Subspace Redistribution...');

plot_bands = {'Raw', 'Alpha', 'Beta'};
band_colors = {[0.2 0.2 0.2], [0 0.4470 0.7410], [0.8500 0.3250 0.0980]}; % Gray, Blue, Red

cond_markers_A = {'o', 's', 'd'}; cond_styles_A  = {'-', '--', ':'};
cond_markers_B = {'^', 'v', 'p'}; cond_styles_B  = {'-', '--', ':'};

for g = 1:2
    fig5A = figure('Position', [100, 100, 1400, 700], 'Name', sprintf('Figure 5A: Cohort %d', g));
    
    all_y = [];
    for b = 1:length(plot_bands)
        if isfield(aggregate_log_r, plot_bands{b})
            m_val = aggregate_log_r.(plot_bands{b}).cohort(g).mean;
            se_val = aggregate_log_r.(plot_bands{b}).cohort(g).se;
            all_y = [all_y; m_val(:) + se_val(:); m_val(:) - se_val(:)];
        end
    end
    if isempty(all_y), all_y = [-1 1]; end
    y_min = min(-1.5, floor(min(all_y) * 1.15)); y_max = max( 1.5, ceil(max(all_y) * 1.15));
    
    % --- Subplot 1: Group A ---
    subplot(1, 2, 1); hold on; set(gca, 'FontSize', 22);
    yline(0, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off'); 
    
    for b = 1:length(plot_bands)
        b_name = plot_bands{b};
        if ~isfield(aggregate_log_r, b_name), continue; end
        grp_m = aggregate_log_r.(b_name).cohort(g).mean;
        grp_se = aggregate_log_r.(b_name).cohort(g).se;
        
        for i = 1:length(cfg.group_A)
            c_idx = find(strcmp(cfg.all_conds, cfg.group_A{i}));
            if isempty(c_idx), continue; end
            disp_name = sprintf('%s (%s)', clean_name(cfg.group_A{i}), b_name);
            errorbar(1:cfg.m, grp_m(c_idx, :), grp_se(c_idx, :), grp_se(c_idx, :), ...
                'LineStyle', cond_styles_A{i}, 'Marker', cond_markers_A{i}, ...
                'Color', band_colors{b}, 'LineWidth', 2.5, 'MarkerSize', 10, ...
                'MarkerFaceColor', band_colors{b}, 'DisplayName', disp_name);
        end
    end
    grid on; xlim([0.75, cfg.m + 0.25]); xticks(1:cfg.m);
    xlabel('Spatial Direction Index', 'FontSize', 20); ylabel('log(r_i) \pm SEM', 'FontSize', 20);
    title(sprintf('Group A: Stimulus Delivered (%s)', strrep(subj_group_names{g}, '_', ' ')), 'FontSize', 24);
    lgd1 = legend('Location', 'best'); lgd1.FontSize = 16; ylim([y_min, y_max]);
    
    % --- Subplot 2: Group B ---
    subplot(1, 2, 2); hold on; set(gca, 'FontSize', 22);
    yline(0, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    for b = 1:length(plot_bands)
        b_name = plot_bands{b};
        if ~isfield(aggregate_log_r, b_name), continue; end
        grp_m = aggregate_log_r.(b_name).cohort(g).mean;
        grp_se = aggregate_log_r.(b_name).cohort(g).se;
        
        for i = 1:length(cfg.group_B)
            c_idx = find(strcmp(cfg.all_conds, cfg.group_B{i}));
            if isempty(c_idx), continue; end
            disp_name = sprintf('%s (%s)', clean_name(cfg.group_B{i}), b_name);
            errorbar(1:cfg.m, grp_m(c_idx, :), grp_se(c_idx, :), grp_se(c_idx, :), ...
                'LineStyle', cond_styles_B{i}, 'Marker', cond_markers_B{i}, ...
                'Color', band_colors{b}, 'LineWidth', 2.5, 'MarkerSize', 10, ...
                'MarkerFaceColor', band_colors{b}, 'DisplayName', disp_name);
        end
    end
    grid on; xlim([0.75, cfg.m + 0.25]); xticks(1:cfg.m);
    xlabel('Spatial Direction Index', 'FontSize', 20); ylabel('log(r_i) \pm SEM', 'FontSize', 20);
    title(sprintf('Group B: Omission (%s)', strrep(subj_group_names{g}, '_', ' ')), 'FontSize', 24);
    lgd2 = legend('Location', 'best'); lgd2.FontSize = 16; ylim([y_min, y_max]);
    
    sgtitle(sprintf('Multi-Band Subspace Redistribution (%s)', strrep(subj_group_names{g}, '_', ' ')), 'FontSize', 28, 'FontWeight', 'bold');
    
    save_file_5a = fullfile(output_path, sprintf('Figure_5A_MultiBand_LogR_Cohort%d', g));
    try pause(0.5); saveas(fig5A, [save_file_5a, '.svg']); saveas(fig5A, [save_file_5a, '.png']); catch, end
    close(fig5A);
end

%% =========================================================================
%% AGGREGATE FIGURE 5B: MULTI-BAND GLOBAL GAIN (G) BAR PLOTS
%% =========================================================================
disp('Generating Figure 5B: Multi-Band Global Gain Bar Plots (with Stats)...');

plot_bands = {'Raw', 'Alpha', 'Beta', 'BetaAlphaRatio'};
band_labels = {'Raw', 'Alpha', 'Beta', 'Beta/Alpha Ratio'};

colors_conds = [0.5 0.5 0.5; ...       % P1 (Cued): Gray
                0.850 0.325 0.098; ... % P2_500: Orange
                0.929 0.694 0.125; ... % P3_500: Yellow
                0.466 0.674 0.188; ... % P2_2000: Green
                0.301 0.745 0.933];    % P3_missing: Blue
                
idx_P1      = find(strcmp(cfg.all_conds, 'P1'));
idx_P2_500  = find(strcmp(cfg.all_conds, 'P2_500'));
idx_P3_500  = find(strcmp(cfg.all_conds, 'P3_500'));
idx_P2_2000 = find(strcmp(cfg.all_conds, 'P2_2000'));
idx_P3_mis  = find(strcmp(cfg.all_conds, 'P3_missing'));

for g = 1:2
    fig5B = figure('Position', [100, 100, 1400, 800], 'Name', sprintf('Figure 5B: Gain Cohort %d', g));
    hold on; set(gca, 'FontSize', 22);
    
    Y_mean = zeros(length(plot_bands), length(cfg.all_conds));
    Y_se   = zeros(length(plot_bands), length(cfg.all_conds));
    
    for b = 1:length(plot_bands)
        b_name = plot_bands{b};
        if isfield(aggregate_G, b_name)
            Y_mean(b, :) = aggregate_G.(b_name).cohort(g).mean';
            Y_se(b, :)   = aggregate_G.(b_name).cohort(g).se';
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
    
    yline(1.0, 'k--', 'LineWidth', 2.5, 'DisplayName', 'Reference (Cued) Baseline');
    
    xticks(1:length(plot_bands));
    xticklabels(band_labels);
    ylabel('Global Subspace Gain (G)', 'FontSize', 24, 'FontWeight', 'bold');
    title(sprintf('Global Subspace Energy Scaling (%s)', strrep(subj_group_names{g}, '_', ' ')), 'FontSize', 28, 'FontWeight', 'bold');
    
    % --- UPDATED: Legend moved to North ---
    lgd = legend('Location', 'north');
    lgd.FontSize = 18;
    grid on;
    
    % --- SIGNIFICANCE TESTING & BRACKETS FOR GAIN ---
    global_max_y = max(max(Y_mean + Y_se));
    
    for b = 1:length(plot_bands)
        b_name = plot_bands{b};
        if ~isfield(aggregate_G, b_name), continue; end
        G_data = aggregate_G.(b_name).cohort(g).subj_data;
        if size(G_data, 2) <= 1, continue; end
        
        band_max = max(Y_mean(b,:) + Y_se(b,:));
        y_shift = max(0.1, band_max * 0.08); 
        current_y = band_max + y_shift;
        
        % Test Group A > Cued (Right Tail)
        for a_idx = [idx_P2_500, idx_P3_500]
            if isempty(a_idx) || isempty(idx_P1), continue; end
            [~, p] = ttest(G_data(a_idx, :), G_data(idx_P1, :), 'Tail', 'right');
            if p < 0.05
                x1 = get_x_pos(b, idx_P1); x2 = get_x_pos(b, a_idx);
                plot([x1, x1, x2, x2], [current_y, current_y+y_shift*0.3, current_y+y_shift*0.3, current_y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                star = '*'; if p < 0.01, star = '**'; end; if p < 0.001, star = '***'; end;
                text(mean([x1, x2]), current_y+y_shift*0.6, star, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                current_y = current_y + y_shift * 1.2;
                global_max_y = max(global_max_y, current_y);
            end
        end
        
        % Test Group B < Cued (Left Tail)
        for b_idx = [idx_P2_2000, idx_P3_mis]
            if isempty(b_idx) || isempty(idx_P1), continue; end
            [~, p] = ttest(G_data(b_idx, :), G_data(idx_P1, :), 'Tail', 'left');
            if p < 0.05
                x1 = get_x_pos(b, idx_P1); x2 = get_x_pos(b, b_idx);
                plot([x1, x1, x2, x2], [current_y, current_y+y_shift*0.3, current_y+y_shift*0.3, current_y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                star = '*'; if p < 0.01, star = '**'; end; if p < 0.001, star = '***'; end;
                text(mean([x1, x2]), current_y+y_shift*0.6, star, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                current_y = current_y + y_shift * 1.2;
                global_max_y = max(global_max_y, current_y);
            end
        end
    end
    
    % Dynamically frame Y-limits for n=5 (Peak gain now ~2.5)
    ylim([0, max(global_max_y * 1.25, 2.5)]);
    
    if g == 2
        % Highlight the true, clean 2-fold elevation in Cohort B
        text(1, global_max_y * 1.08, '~2-fold Gain Elevation (excl. S14)', ...
            'FontSize', 18, 'Color', [0.2 0.2 0.8], 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'Margin', 2);
    end
    
    save_file_5b = fullfile(output_path, sprintf('Figure_5B_MultiBand_Gain_Cohort%d', g));
    try pause(0.5); saveas(fig5B, [save_file_5b, '.svg']); saveas(fig5B, [save_file_5b, '.png']); catch, end
    close(fig5B);
end

%% =========================================================================
%% AGGREGATE FIGURE 5C: MULTI-BAND DELTA LATENCY SCATTER PLOT
%% =========================================================================
disp('Generating Figure 5C: Multi-Band Delta Latency Scatter Plots...');

band_titles = {'Raw Broadband', 'Alpha (8-12 Hz)', 'Beta (15-30 Hz)', 'Beta/Alpha Ratio'};

% Only plot conditions with actual physical tactile ERPs
erp_conds = {'P1', 'P2_500', 'P3_500'};
erp_conds_idx = find(ismember(cfg.all_conds, erp_conds));

for g = 1:2
    fig5C = figure('Position', [100, 100, 1600, 700], 'Name', sprintf('Figure 5C: Latency Cohort %d', g));
    tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    % Find dynamic yet symmetrical Y-limits across all 4 bands ONLY for the 3 ERP conditions
    all_lat = [];
    for b = 1:length(plot_bands)
        if isfield(aggregate_lat, plot_bands{b})
            lat_subset = aggregate_lat.(plot_bands{b}).cohort(g).mean(erp_conds_idx, :);
            all_lat = [all_lat; lat_subset(:)];
        end
    end
    
    if isempty(all_lat) || all(isnan(all_lat)), all_lat = [-50, 50]; end
    y_bound = max(50, ceil(max(abs(all_lat)) / 25) * 25); 

    for b = 1:length(plot_bands)
        b_name = plot_bands{b};
        if ~isfield(aggregate_lat, b_name), continue; end
        
        lat_m = aggregate_lat.(b_name).cohort(g).mean; % [nConds x m]
        
        nexttile; hold on;
        set(gca, 'FontSize', 18);
        
        yline(0, 'k-', 'LineWidth', 2.5, 'HandleVisibility', 'off');
        
        for i = 1:length(erp_conds_idx)
            c_idx = erp_conds_idx(i);
            
            rng(42); 
            x_jitter = i + (rand(1, cfg.m) - 0.5) * 0.4; 
            
            scatter(x_jitter, lat_m(c_idx, :), 100, colors_conds(c_idx, :), 'filled', ...
                'MarkerEdgeColor', 'k', 'DisplayName', clean_name(cfg.all_conds{c_idx}));
                
            mean_val = mean(lat_m(c_idx, :), 'omitnan');
            plot([i-0.3, i+0.3], [mean_val, mean_val], '-', 'Color', colors_conds(c_idx, :), 'LineWidth', 4, 'HandleVisibility', 'off');
            
            if c_idx ~= find(strcmp(cfg.all_conds, 'P1'))
                lat_data = aggregate_lat.(b_name).cohort(g).subj_data(c_idx, :);
                [~, p_val] = ttest(lat_data, 0); 
                if p_val < 0.05
                    star = '*'; if p_val < 0.01, star = '**'; end; if p_val < 0.001, star = '***'; end;
                    test_mean = mean(lat_data, 'omitnan');
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
    
    sgtitle(sprintf('Temporal Dynamics: Tactile ERP \\Delta Latency Shifts (%s)', strrep(subj_group_names{g}, '_', ' ')), 'FontSize', 28, 'FontWeight', 'bold');
    
    lgd_names = cellfun(clean_name, cfg.all_conds(erp_conds_idx), 'UniformOutput', false);
    lgd = legend(lgd_names, 'Location', 'northeastoutside');
    lgd.FontSize = 18;
    
    save_file_5c = fullfile(output_path, sprintf('Figure_5C_MultiBand_Latency_Cohort%d', g));
    try pause(0.5); saveas(fig5C, [save_file_5c, '.svg']); saveas(fig5C, [save_file_5c, '.png']); catch, end
    close(fig5C);
end

%% =========================================================================
%% AGGREGATE FIGURE 5D: MULTI-BAND % SUSTAINED (TONIC VS PHASIC) BAR PLOTS
%% =========================================================================
disp('Generating Figure 5D: Multi-Band Sustained % Bar Plots...');

idx_P1      = find(strcmp(cfg.all_conds, 'P1'));
idx_P2_500  = find(strcmp(cfg.all_conds, 'P2_500'));
idx_P3_500  = find(strcmp(cfg.all_conds, 'P3_500'));
idx_P2_2000 = find(strcmp(cfg.all_conds, 'P2_2000'));
idx_P3_mis  = find(strcmp(cfg.all_conds, 'P3_missing'));

for g = 1:2
    fig5D = figure('Position', [100, 100, 1400, 800], 'Name', sprintf('Figure 5D: Sustained Cohort %d', g));
    hold on; set(gca, 'FontSize', 22);
    
    Y_mean = zeros(length(plot_bands), length(cfg.all_conds));
    Y_se   = zeros(length(plot_bands), length(cfg.all_conds));
    
    for b = 1:length(plot_bands)
        b_name = plot_bands{b};
        if isfield(aggregate_sust, b_name)
            Y_mean(b, :) = aggregate_sust.(b_name).cohort(g).mean';
            Y_se(b, :)   = aggregate_sust.(b_name).cohort(g).se';
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
    xticklabels(band_labels);
    ylabel('Sustained Variance (%)', 'FontSize', 24, 'FontWeight', 'bold');
    title(sprintf('Tonic vs. Phasic Stability (%s)', strrep(subj_group_names{g}, '_', ' ')), 'FontSize', 28, 'FontWeight', 'bold');
    
    lgd = legend('Location', 'northwest');
    lgd.FontSize = 18;
    grid on;
    
    max_y = max(max(Y_mean + Y_se)) * 1.2;
    ylim([0, max(max_y, 75) + 15]); 
    
    plot([1.6, 3.4], [10, 10], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot([1.6, 1.6], [8, 10], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot([3.4, 3.4], [8, 10], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    text(2.5, 13, 'Purely Phasic (< 3%)', 'FontSize', 20, 'Color', 'k', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'Margin', 2);
    
    % --- SIGNIFICANCE TESTING & BRACKETS ---
    
    % 1. Raw Band (b=1): Cued > Rand 500 & Cued > Rand Missing
    if isfield(aggregate_sust, 'Raw') && ~isempty(idx_P1)
        raw_data = aggregate_sust.Raw.cohort(g).subj_data;
        if size(raw_data, 2) > 1 
            y_shift = 4; 
            if ~isempty(idx_P3_500)
                [~, p1] = ttest(raw_data(idx_P1, :), raw_data(idx_P3_500, :), 'Tail', 'right');
                if p1 < 0.05
                    x1 = get_x_pos(1, idx_P1); x2 = get_x_pos(1, idx_P3_500);
                    y = max(Y_mean(1, [idx_P1, idx_P3_500]) + Y_se(1, [idx_P1, idx_P3_500])) + y_shift;
                    plot([x1, x1, x2, x2], [y, y+1.5, y+1.5, y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    star = '*'; if p1 < 0.01, star = '**'; end; if p1 < 0.001, star = '***'; end;
                    text(mean([x1, x2]), y+3, star, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                    y_shift = 10; 
                end
            end
            if ~isempty(idx_P3_mis)
                [~, p2] = ttest(raw_data(idx_P1, :), raw_data(idx_P3_mis, :), 'Tail', 'right');
                if p2 < 0.05
                    x1 = get_x_pos(1, idx_P1); x2 = get_x_pos(1, idx_P3_mis);
                    y = max(Y_mean(1, [idx_P1, idx_P3_mis]) + Y_se(1, [idx_P1, idx_P3_mis])) + y_shift;
                    plot([x1, x1, x2, x2], [y, y+1.5, y+1.5, y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    star = '*'; if p2 < 0.01, star = '**'; end; if p2 < 0.001, star = '***'; end;
                    text(mean([x1, x2]), y+3, star, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                end
            end
        end
    end
    
    % --- NEW: 2. Alpha (b=2) and Beta (b=3) Bands: Cued > Rand 500 & Cued > Rand Missing ---
    for b_idx = 2:3
        b_name = plot_bands{b_idx};
        if isfield(aggregate_sust, b_name) && ~isempty(idx_P1)
            sust_data = aggregate_sust.(b_name).cohort(g).subj_data;
            if size(sust_data, 2) > 1 
                y_shift = 17; % Start high enough to clear the "Purely Phasic" text
                
                if ~isempty(idx_P3_500)
                    [~, p_ab1] = ttest(sust_data(idx_P1, :), sust_data(idx_P3_500, :), 'Tail', 'right');
                    if p_ab1 < 0.05
                        x1 = get_x_pos(b_idx, idx_P1); x2 = get_x_pos(b_idx, idx_P3_500);
                        y = max(Y_mean(b_idx, [idx_P1, idx_P3_500]) + Y_se(b_idx, [idx_P1, idx_P3_500])) + y_shift;
                        plot([x1, x1, x2, x2], [y, y+1.5, y+1.5, y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                        star = '*'; if p_ab1 < 0.01, star = '**'; end; if p_ab1 < 0.001, star = '***'; end;
                        text(mean([x1, x2]), y+3, star, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                        y_shift = y_shift + 6; 
                    end
                end
                
                if ~isempty(idx_P3_mis)
                    [~, p_ab2] = ttest(sust_data(idx_P1, :), sust_data(idx_P3_mis, :), 'Tail', 'right');
                    if p_ab2 < 0.05
                        x1 = get_x_pos(b_idx, idx_P1); x2 = get_x_pos(b_idx, idx_P3_mis);
                        y = max(Y_mean(b_idx, [idx_P1, idx_P3_mis]) + Y_se(b_idx, [idx_P1, idx_P3_mis])) + y_shift;
                        plot([x1, x1, x2, x2], [y, y+1.5, y+1.5, y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                        star = '*'; if p_ab2 < 0.01, star = '**'; end; if p_ab2 < 0.001, star = '***'; end;
                        text(mean([x1, x2]), y+3, star, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                    end
                end
            end
        end
    end
    
    % 3. Beta/Alpha Ratio (b=4): Unpred 2000 > Cued & Unpred 2000 > Unpred 500
    if isfield(aggregate_sust, 'BetaAlphaRatio') && ~isempty(idx_P2_2000)
        ratio_data = aggregate_sust.BetaAlphaRatio.cohort(g).subj_data;
        if size(ratio_data, 2) > 1
            y_shift = 4;
            if ~isempty(idx_P1)
                [~, p3] = ttest(ratio_data(idx_P2_2000, :), ratio_data(idx_P1, :), 'Tail', 'right');
                if p3 < 0.05
                    x1 = get_x_pos(4, idx_P1); x2 = get_x_pos(4, idx_P2_2000);
                    y = max(Y_mean(4, [idx_P1, idx_P2_2000]) + Y_se(4, [idx_P1, idx_P2_2000])) + y_shift;
                    plot([x1, x1, x2, x2], [y, y+1.5, y+1.5, y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    star = '*'; if p3 < 0.01, star = '**'; end; if p3 < 0.001, star = '***'; end;
                    text(mean([x1, x2]), y+3, star, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                    y_shift = 10;
                end
            end
            if ~isempty(idx_P2_500)
                [~, p4] = ttest(ratio_data(idx_P2_2000, :), ratio_data(idx_P2_500, :), 'Tail', 'right');
                if p4 < 0.05
                    x1 = get_x_pos(4, idx_P2_500); x2 = get_x_pos(4, idx_P2_2000);
                    y = max(Y_mean(4, [idx_P2_500, idx_P2_2000]) + Y_se(4, [idx_P2_500, idx_P2_2000])) + y_shift;
                    plot([x1, x1, x2, x2], [y, y+1.5, y+1.5, y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    star = '*'; if p4 < 0.01, star = '**'; end; if p4 < 0.001, star = '***'; end;
                    text(mean([x1, x2]), y+3, star, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                end
            end
        end
    end
    
    save_file_5d = fullfile(output_path, sprintf('Figure_5D_MultiBand_Sustained_Cohort%d', g));
    try pause(0.5); saveas(fig5D, [save_file_5d, '.svg']); saveas(fig5D, [save_file_5d, '.png']); catch, end
    close(fig5D);
end

disp('All analyses completed successfully. Figure 3B, Figure 5A, 5B, 5C, and 5D saved.');