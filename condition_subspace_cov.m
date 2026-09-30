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

% Loop through the conditions
for c = 1:length(conditions)
    condition = conditions{c};
    in_dir = fullfile(input_path, condition);
    set_files = dir(fullfile(in_dir, '*.set'));
    names_sorted = sort(cellstr({set_files.name})');
    
    if isempty(names_sorted)
        fprintf('No .set files for %s — skipping\n', condition);
        continue;
    end
    
    if strcmp(condition, 'P2')
        excel_file_path = fullfile(input_path, 'Indexes for P2.xlsx');
        epoch_trials_p2_500ms = readmatrix(excel_file_path, 'Sheet', 'Audio onset with 500 ms tactile');
        epoch_trials_p2_2000ms = readmatrix(excel_file_path, 'Sheet', 'Audio onset with 2000 ms tactil');
    elseif strcmp(condition, 'P3')
        excel_file_path = fullfile(input_path, 'Indexes for P3.xlsx');
        epoch_trials_p3_500ms = readmatrix(excel_file_path, 'Sheet', 'Audio onset with 500 ms tactile');
        epoch_trials_p3_missing = readmatrix(excel_file_path, 'Sheet', 'Audio onset with missing tactil');
    end
    
    for s = 1:length(names_sorted)
        file_to_load = names_sorted{s}; 
        fprintf('Loading Subj %d/%d (%s) for condition: %s\n', s, length(names_sorted), file_to_load, condition);
        
        EEG = pop_loadset('filename', file_to_load, 'filepath', in_dir);
        if isempty(fs), fs = EEG.srate; end
        
        if isempty(time_ms_eeg)
            time_ms_eeg = linspace(0, 3.5 , size(EEG.data, 2)); 
        end
        
        if ismember(condition, {'BLA', 'BLT', 'P1'})
            epoch_trials = 1:2:EEG.trials;
            data_all_conds.Raw.(condition){s} = EEG.data(:, :, epoch_trials);
        elseif strcmp(condition, 'P2')
            epoch_trials_p2 = 1:2:EEG.trials;     
            data_all_conds.Raw.(condition){s}   = EEG.data(:, :, epoch_trials_p2);
            data_all_conds.Raw.P2_500{s}        = EEG.data(:, :, epoch_trials_p2_500ms);
            data_all_conds.Raw.P2_2000{s}       = EEG.data(:, :, epoch_trials_p2_2000ms);
        elseif strcmp(condition, 'P3')
            epoch_trials_p3 = 1:2:EEG.trials;
            data_all_conds.Raw.(condition){s}   = EEG.data(:, :, epoch_trials_p3);
            data_all_conds.Raw.P3_500{s}        = EEG.data(:, :, epoch_trials_p3_500ms);
            data_all_conds.Raw.P3_missing{s}    = EEG.data(:, :, epoch_trials_p3_missing);
        end
    end
end

%% 0.2 Apply Alpha and Beta Zero-Phase Filtering & Compute Beta/Alpha Ratio
disp('Applying Zero-Phase FIR Filters for Alpha and Beta...');
fn = fs / 2; 
ord_alpha = round(0.250 * fs); % 250 ms
b_alpha   = fir1(ord_alpha, [8 12] / fn, 'bandpass');

ord_beta  = round(0.125 * fs); % 125 ms
b_beta    = fir1(ord_beta, [13 30] / fn, 'bandpass');

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
if max(time_ms_eeg) > 2.0
    time_s = time_ms_eeg - 1.000; 
else
    time_s = time_ms_eeg;
end
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

% --- NEW: Aggregate Struct for Figure 3B Scree Plot ---
aggregate_spectra = struct();

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
    
    % --- NEW: Capture all 32 components for cumulative plot ---
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
        
        % --- NEW: Save full spectrum fraction ---
        lam_full_all(:, s)    = lam / sum(lam); 
    end
    
    % --- NEW: Save aggregate cumulative variance (mean across subjects) in percentage ---
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
    % Group arrays: [Conditions x Metric/Dims x Reps x Subjects]
    G_splits     = nan(length(cfg.all_conds), cfg.nRep, num_subjects);
    log_r_splits = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
    lat_splits   = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
    sust_splits  = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
    traj_splits  = nan(length(cfg.all_conds), cfg.m, num_windows, cfg.nRep, num_subjects);

    for subj = 1:num_subjects
        subj_dir = fullfile(output_path, sprintf('subj%02d', subj), band_name, cfg.ref_cond);
        if ~exist(subj_dir, 'dir'), mkdir(subj_dir); end
        save(fullfile(subj_dir, 'cfg.mat'), 'cfg');
        
        % Stable alignment basis for this subject
        Xi_stable = Xi_aligned(:, :, subj);
        
        for rep = 1:cfg.nRep
            rng(cfg.seed + rep + subj*1000); 
            
            % 1. Dynamic Split-Half
            raw_ref = current_data.(cfg.ref_cond){subj};
            nRef = size(raw_ref, 3);
            perm = randperm(nRef);
            iA = perm(1:cfg.nTrialMatch);
            iB = perm(cfg.nTrialMatch + (1:cfg.nTrialMatch));
            
            % Half A Basis recalculation per rep
            erpA = mean(raw_ref(:, :, iA), 3, 'omitnan');
            bA   = mean(erpA(:, iBase), 2);
            
            % Calculate unregularized centered covariance
            C_A_cent = covCentered(erpA(:, iWin), bA);
            
            [V, D] = eig(C_A_cent, 'vector');
            [lam, ord] = sort(D, 'descend');
            V = V(:, ord);
            Xi_rep = V(:, 1:cfg.m);
            lam_frac_rep = lam(1:cfg.m) / sum(lam);
            gap_rep = (lam(1:cfg.m-1) - lam(2:cfg.m)) ./ lam(1:cfg.m-1);
            
            % Align dynamic rep basis to stable subject basis
            for k = 1:cfg.m
                if dot(Xi_rep(:, k), Xi_stable(:, k)) < 0
                    Xi_rep(:, k) = -Xi_rep(:, k);
                end
            end
            
            basis_vars = struct('Xi', Xi_rep, 'lam', lam, 'lam_frac', lam_frac_rep, 'gap', gap_rep);
            
            % Half B Ceiling Evaluation
            erpB = mean(raw_ref(:, :, iB), 3, 'omitnan');
            bB   = mean(erpB(:, iBase), 2);
            XB   = erpB(:, iWin) - bB;
            
            v_ceil = sum((Xi_rep' * XB).^2, 2) / T_win;
            V_ceil = sum(v_ceil);
            
            % Calculate Ceiling Latency for ALL m components
            proj_B = Xi_rep' * XB; % [m x T]
            [~, max_idx_B] = max(abs(proj_B), [], 2);
            t_lat_ceil = t_eval(max_idx_B); % [m x 1]
            
            % 2. Projections
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
                % Static 0.0-0.5s Projection
                v_k = sum((Xi_rep' * Xk).^2, 2) / T_win;
                V_k = sum(v_k);
                G_k = V_k / V_ceil;
                r_k = v_k ./ v_ceil;
                log_r_k = log(r_k);
                
                proj_vars.(cond) = struct('v_k', v_k, 'G', G_k, 's_k', v_k/V_k, 'r', r_k, 'log_r', log_r_k);
                
                % Store in Group Arrays
                G_splits(c, rep, subj)        = G_k;
                log_r_splits(c, :, rep, subj) = log_r_k;
                
                % Calculate Projected Latency for ALL m components
                proj_k = Xi_rep' * Xk; % [m x T]
                [~, max_idx_k] = max(abs(proj_k), [], 2);
                lat_splits(c, :, rep, subj) = (t_eval(max_idx_k) - t_lat_ceil)' * 1000;
                
                xbar = mean(Xk, 2); Xd = Xk - xbar;
                v_sust = (Xi_rep' * xbar).^2;
                v_dyn  = sum((Xi_rep' * Xd).^2, 2) / T_win;
                sust_splits(c, :, rep, subj) = v_sust ./ (v_sust + v_dyn + eps);

                % --- Moving Window Dynamic Projection ---
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
    
    % 1. Calculate TRUE empirical mean pairwise absolute cosine similarity
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
    
    % 2. Generate NULL distributions via Channel Shuffling
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
    
    % 3. Calculate p-values and 95% Confidence Thresholds
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
    
    % 4. Plot Null Distributions vs True Values
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
    subj_mean_log_r = squeeze(mean(log_r_splits, 3, 'omitnan')); % [nCond x m x Subjs]
    subj_mean_traj  = squeeze(mean(traj_splits, 4, 'omitnan'));  % [nCond x m x nWin x Subjs]
    
    % 2. Define the two subject subgroups
    g1_subjs = intersect(4:11, 1:num_subjects); % Protects against out-of-bounds
    g2_subjs = setdiff(1:num_subjects, g1_subjs);
    
    subj_groups = {g1_subjs, g2_subjs};
    subj_group_names = {'Subjs_4_to_11', 'Other_Subjs'};
    
    % =========================================================================
    % 3. PRE-CALCULATE SYNCHRONIZED Y-LIMITS & C-LIMITS ACROSS BOTH COHORTS
    % =========================================================================
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
            
            % Evaluate Cohort 1
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
            
            % Evaluate Cohort 2
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
            
            % --- UPDATED: Reduced Checkpoint 3a legend font size to 20 ---
            if k_idx == 1, lgdA = legend('Location', 'best'); lgdA.FontSize = 20; end
            
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
            
            % --- UPDATED: Reduced Checkpoint 3a legend font size to 20 ---
            if k_idx == 1, lgdB = legend('Location', 'best'); lgdB.FontSize = 20; end
            
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

% Truncation highlight at m=6
xline(cfg.m, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');

% Highlight >85% at m=6
yline(85, 'r:', 'LineWidth', 2, 'DisplayName', '85% Variance Threshold');
plot(cfg.m, 85, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');

% Annotation for Beta/Alpha Ratio at m=2
xline(2, ':', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2, 'HandleVisibility', 'off');
text(2.2, 50, 'Beta/Alpha intrinsic 2D manifold (p < 0.05)', ...
    'FontSize', 18, 'Color', [0.4940 0.1840 0.5560], 'FontWeight', 'bold', 'BackgroundColor', 'w');
    
text(cfg.m + 0.2, 70, sprintf('Truncation at m=%d\n(>85%% variance)', cfg.m), ...
    'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold', 'BackgroundColor', 'w');
    
xlim([0.5, 15]); % Limit x-axis to 15 to show the elbow clearly
xticks(1:15);
ylim([0, 100]);
yticks(0:10:100);

xlabel('Number of Components (Spatial Directions)', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Cumulative Variance Explained (%)', 'FontSize', 24, 'FontWeight', 'bold');
title('Intrinsic Dimensionality Across Frequency Bands', 'FontSize', 28, 'FontWeight', 'bold');

lgd = legend('Location', 'southeast');
lgd.FontSize = 20;
grid on;

% --- EXPORTING AS SVG AND PNG ---
save_file_3b = fullfile(output_path, 'Figure_3B_Aggregate_ScreePlot');
try 
    pause(0.5); 
    saveas(fig3B, [save_file_3b, '.svg']); 
    saveas(fig3B, [save_file_3b, '.png']); 
catch
    warning('Could not save %s.', save_file_3b); 
end
close(fig3B);

disp('All analyses completed successfully. Figure 3B saved.');