clear; clc;
% --- paths ---
if exist('I:\', 'dir')
    input_path = 'I:\My Drive\Data\New Data\EEG epoched\';
    base_output_path = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper\ch-ch_cov';
elseif exist('H:\', 'dir')
    input_path = 'H:\My Drive\Data\New Data\EEG epoched\';
    base_output_path = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper\ch-ch_cov';
elseif exist('G:\', 'dir')
    input_path = 'G:\My Drive\Data\New Data\EEG epoched\';
    base_output_path = 'C:\Users\sdabiri\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper\ch-ch_cov';
else
    error('Unknown system: Cannot determine input and output paths.');
end
output_path = fullfile(base_output_path);

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
    
    % Load excel file paths for P2 and P3 ONCE per condition
    if strcmp(condition, 'P2')
        excel_file_path = fullfile(input_path, 'Indexes for P2.xlsx');
        epoch_trials_p2_500ms = readmatrix(excel_file_path, 'Sheet', 'Audio onset with 500 ms tactile');
        epoch_trials_p2_2000ms = readmatrix(excel_file_path, 'Sheet', 'Audio onset with 2000 ms tactil');
    elseif strcmp(condition, 'P3')
        excel_file_path = fullfile(input_path, 'Indexes for P3.xlsx');
        epoch_trials_p3_500ms = readmatrix(excel_file_path, 'Sheet', 'Audio onset with 500 ms tactile');
        epoch_trials_p3_missing = readmatrix(excel_file_path, 'Sheet', 'Audio onset with missing tactil');
    end

    % Load ALL subjects for this condition
    for s = 1:length(names_sorted)
        file_to_load = names_sorted{s}; 
        fprintf('Loading Subj %d/%d (%s) for condition: %s\n', s, length(names_sorted), file_to_load, condition);
        
        EEG = pop_loadset('filename', file_to_load, 'filepath', in_dir);
        if isempty(fs), fs = EEG.srate; end
        
        if isempty(time_ms_eeg)
            time_ms_eeg = linspace(0, 3.5 , size(EEG.data, 2)); 
        end
        
        % Extract data and store in Raw struct
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

%% 0.2 Apply Alpha and Beta Zero-Phase Filtering
disp('Applying FIR Filters for Alpha (8-12 Hz) and Beta (13-30 Hz)...');

% Design FIR Filters
fn = fs / 2; % Nyquist frequency
ord_alpha = round(0.250 * fs); % 250 ms window
b_alpha   = fir1(ord_alpha, [8 12] / fn, 'bandpass');

ord_beta  = round(0.125 * fs); % 125 ms window
b_beta    = fir1(ord_beta, [13 30] / fn, 'bandpass');

cond_names = fieldnames(data_all_conds.Raw);
for c = 1:length(cond_names)
    c_name = cond_names{c};
    for s = 1:length(data_all_conds.Raw.(c_name))
        raw_data = data_all_conds.Raw.(c_name){s};
        if isempty(raw_data), continue; end
        
        [nc, nt, ntr] = size(raw_data);
        
        % Permute to [Time x Channels x Trials] and reshape to 2D [Time x (Channels*Trials)]
        % This ensures filtfilt processes each channel of each trial independently
        X_perm = permute(raw_data, [2, 1, 3]);
        X_2D   = reshape(X_perm, nt, nc * ntr);
        
        % Filter Alpha & Reshape Back
        X_alpha_2D = filtfilt(b_alpha, 1, double(X_2D));
        data_all_conds.Alpha.(c_name){s} = permute(reshape(X_alpha_2D, nt, nc, ntr), [2, 1, 3]);
        
        % Filter Beta & Reshape Back
        X_beta_2D  = filtfilt(b_beta, 1, double(X_2D));
        data_all_conds.Beta.(c_name){s}  = permute(reshape(X_beta_2D, nt, nc, ntr), [2, 1, 3]);
    end
end
disp('Filtering complete.');

%% 1. Configuration Parameters
cfg = struct();
cfg.num_ch        = 32;
cfg.fs            = fs;                
cfg.tStart        = -1.000;            
cfg.win           = [0.000, 0.500];    
cfg.base          = [-1.000, -0.800];  
cfg.base_robust   = [-0.800, -0.600];  
cfg.nTrialMatch   = 30;                
cfg.nRep          = 50;                
cfg.m             = 6;                 
cfg.alpha         = 0.05;              
cfg.seed          = 20260910;

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

%% 3. Helper Functions
shrinkCov = @(C, a) (1 - a) * C + a * (trace(C) / size(C, 1)) * eye(size(C, 1));
covCentered = @(X, b, a) shrinkCov((( (X - b) - mean(X - b, 2) ) * ( (X - b) - mean(X - b, 2) )') / size(X, 2), a);


%% =========================================================================
%% MASTER LOOP OVER FREQUENCY BANDS
%% =========================================================================
bands_to_process = {'Raw', 'Alpha', 'Beta'};

for band_idx = 1:length(bands_to_process)
    band_name = bands_to_process{band_idx};
    current_data = data_all_conds.(band_name); % Extract specific band data
    
    fprintf('\n\n======================================================\n');
    fprintf('       PROCESSING SIGNAL BAND: %s\n', upper(band_name));
    fprintf('======================================================\n');

    % ----------------------------------------------------------------------
    %% 4. CHECKPOINT 1: Multi-Subject Eigenbasis Extraction (BLT Ref)
    % ----------------------------------------------------------------------
    num_subjects = length(current_data.BLT);
    Xi_all_subjs = nan(cfg.num_ch, cfg.m, num_subjects);
    lam_frac_all = nan(cfg.m, num_subjects);
    gap_all      = nan(cfg.m - 1, num_subjects);
    
    for s = 1:num_subjects
        raw_BLT = current_data.BLT{s};
        n_trials = size(raw_BLT, 3);
        if n_trials < (cfg.nTrialMatch * 2)
            error('Subject %d has insufficient BLT trials (%d < 60)', s, n_trials);
        end
        
        rng(cfg.seed + s);
        perm = randperm(n_trials);
        iA = perm(1:cfg.nTrialMatch);
        
        erpA = mean(raw_BLT(:, :, iA), 3, 'omitnan');
        bA   = mean(erpA(:, iBase), 2);
        
        C_A_cent = covCentered(erpA(:, iWin), bA, cfg.alpha);
        
        [V, D] = eig(C_A_cent, 'vector');
        [lam, ord] = sort(D, 'descend');
        V = V(:, ord);
        
        Xi_all_subjs(:, :, s) = V(:, 1:cfg.m);
        lam_frac_all(:, s)   = lam(1:cfg.m) / sum(lam);
        gap_all(:, s)        = (lam(1:cfg.m-1) - lam(2:cfg.m)) ./ lam(1:cfg.m-1);
    end
    
    Xi_aligned = Xi_all_subjs;
    for s = 2:num_subjects
        for k = 1:cfg.m
            if dot(Xi_aligned(:, k, s), Xi_aligned(:, k, 1)) < 0
                Xi_aligned(:, k, s) = -Xi_aligned(:, k, s);
            end
        end
    end
    Xi_bar = mean(Xi_aligned, 3);
    
    % --- Deliverable A ---
    figA = figure('Position', [100, 100, 1000, 400], 'Name', sprintf('Chk 1 (BLT) - A (%s)', band_name));
    subplot(1, 2, 1);
    plot(1:cfg.m, lam_frac_all, '-o', 'LineWidth', 1.5);
    grid on; xlabel('Direction Index', 'FontSize', 16); ylabel('Fraction of Total Variance', 'FontSize', 16);
    title(sprintf('Deliverable A: Eigenvalue Spectra [%s]', band_name), 'FontSize', 16);
    
    subplot(1, 2, 2);
    bar(mean(gap_all, 2)); hold on;
    errorbar(1:cfg.m-1, mean(gap_all, 2), std(gap_all, 0, 2), 'k.', 'LineWidth', 1.2);
    grid on; xlabel('Gap Index', 'FontSize', 16); ylabel('Relative Gap', 'FontSize', 16);
    title(sprintf('Deliverable A: Mean Subspace Gaps [%s]', band_name), 'FontSize', 16);
    
    % --- Deliverable B ---
    n_rows = min(num_subjects, 4);
    fig_w = 200 * cfg.m + 140; 
    fig_h = 190 * n_rows + 40; 
    figure('Position', [50, 50, fig_w, fig_h], 'Name', sprintf('Chk 1 (BLT) - B (%s)', band_name));
    for s = 1:n_rows
        for k = 1:cfg.m
            ax = subplot(n_rows, cfg.m, (s-1)*cfg.m + k);
            pos = get(ax, 'Position');
            set(ax, 'Position', [pos(1) * 0.90 + 0.03, pos(2) * 0.95, pos(3) * 0.90, pos(4) * 0.95]);
            topoplot(Xi_aligned(:, k, s), EEG.chanlocs, 'numcontour', 0);
            clim([-0.4 0.4]); colormap('jet');
            if s == 1
                text(0.5, 1.18, sprintf('\\xi_%d', k), 'Units', 'normalized', 'FontSize', 22, ...
                    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            end
            if k == 1
                text(-0.28, 0.5, sprintf('Subj %d', s), 'Units', 'normalized', 'FontSize', 18, ...
                    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
            end
        end
    end
    cb = colorbar('Position', [0.88, 0.12, 0.022, 0.72]); 
    cb.FontSize = 13; cb.FontWeight = 'bold'; cb.Ticks = [-0.4, -0.2, 0, 0.2, 0.4];
    cb.Label.String = 'Spatial Weight (a.u.)'; cb.Label.FontSize = 15; cb.Label.FontWeight = 'bold';
    
    % --- Deliverable C ---
    sim_matrix = nan(num_subjects, cfg.m);
    for s = 1:num_subjects
        for k = 1:cfg.m
            sim_matrix(s, k) = abs(dot(Xi_aligned(:, k, s), Xi_bar(:, k))) / (norm(Xi_aligned(:, k, s)) * norm(Xi_bar(:, k)));
        end
    end
    figure('Position', [150, 150, 600, 400], 'Name', sprintf('Chk 1 (BLT) - C (%s)', band_name));
    bar(mean(sim_matrix, 1)); hold on;
    errorbar(1:cfg.m, mean(sim_matrix, 1), std(sim_matrix, 0, 1), 'k.', 'LineWidth', 1.5);
    grid on; xlabel('Spatial Direction Index', 'FontSize', 16); ylabel('Cosine Similarity', 'FontSize', 16);
    title(sprintf('Deliverable C: Consistency Across Subjects [%s]', band_name), 'FontSize', 16); ylim([0 1]);

    % ----------------------------------------------------------------------
    %% 5. PROJECTIONS & METRICS (BLT Ref)
    % ----------------------------------------------------------------------
    cfg.ref_cond  = 'BLT';
    cfg.group_A   = {'P1', 'P2_500', 'P3_500'};        
    cfg.group_B   = {'P2_2000', 'P3_missing'};         
    cfg.all_conds = [{'BLT'}, cfg.group_A, cfg.group_B];

    target_subj = 1; 
    Xi = Xi_aligned(:, :, target_subj);
    
    G_splits       = nan(length(cfg.all_conds), cfg.nRep);
    log_r_splits   = nan(length(cfg.all_conds), cfg.m, cfg.nRep);
    lat_splits     = nan(length(cfg.all_conds), cfg.nRep);
    sust_splits    = nan(length(cfg.all_conds), cfg.m, cfg.nRep);
    
    for rep = 1:cfg.nRep
        rng(cfg.seed + rep);
        nB = size(current_data.BLT{target_subj}, 3);
        perm = randperm(nB);
        iB = perm(cfg.nTrialMatch + (1:cfg.nTrialMatch));
        
        erpB = mean(current_data.BLT{target_subj}(:, :, iB), 3, 'omitnan');
        bB   = mean(erpB(:, iBase), 2);
        XB   = erpB(:, iWin) - bB;
        
        v_ceil = sum((Xi' * XB).^2, 2) / T_win;
        V_ceil = sum(v_ceil);
        proj_B = Xi(:, 1)' * XB;
        [~, max_idx_B] = max(abs(proj_B));
        t_lat_ceil = t_eval(max_idx_B);
        
        for c = 1:length(cfg.all_conds)
            cond = cfg.all_conds{c};
            if ~isfield(current_data, cond), continue; end
            
            raw_cond = current_data.(cond){target_subj};
            if strcmp(cond, 'BLT')
                Xk = XB; bk = bB;
            else
                idx_k = randperm(size(raw_cond, 3), cfg.nTrialMatch);
                erp_k = mean(raw_cond(:, :, idx_k), 3, 'omitnan');
                bk    = mean(erp_k(:, iBase), 2);
                Xk    = erp_k(:, iWin) - bk;
            end
            
            v_k = sum((Xi' * Xk).^2, 2) / T_win;
            V_k = sum(v_k);
            G_splits(c, rep)        = V_k / V_ceil;
            log_r_splits(c, :, rep) = log(v_k ./ v_ceil);
            
            proj_k = Xi(:, 1)' * Xk;
            [~, max_idx_k] = max(abs(proj_k));
            lat_splits(c, rep) = (t_eval(max_idx_k) - t_lat_ceil) * 1000;
            
            xbar = mean(Xk, 2); Xd = Xk - xbar;
            v_sust = (Xi' * xbar).^2;
            v_dyn  = sum((Xi' * Xd).^2, 2) / T_win;
            sust_splits(c, :, rep) = v_sust ./ (v_sust + v_dyn + eps);
        end
    end

    % ----------------------------------------------------------------------
    %% 6. CHECKPOINT 2: Results Display (BLT Ref)
    % ----------------------------------------------------------------------
    mean_log_r = mean(log_r_splits, 3, 'omitnan'); 
    prc_log_r  = prctile(log_r_splits, [5 95], 3); 
    
    all_low  = prc_log_r(:, :, 1);
    all_high = prc_log_r(:, :, 2);
    y_min    = min(-1.5, floor(min(all_low(:)) * 1.15));
    y_max    = max( 1.5, ceil(max(all_high(:)) * 1.15));
    y_limits = [y_min, y_max]; 
    
    figure('Position', [100, 100, 1200, 500], 'Name', sprintf('Chk 2 (BLT) - %s', band_name));
    
    subplot(1, 2, 1); hold on;
    yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ceiling (BLT B)');
    colors_A = {[0.85 0.32 0.09], [0.92 0.69 0.12], [0.49 0.18 0.55]};
    for i = 1:length(cfg.group_A)
        c_idx = find(strcmp(cfg.all_conds, cfg.group_A{i}));
        if isempty(c_idx), continue; end
        y_val = mean_log_r(c_idx, :);              
        err_low = y_val - prc_log_r(c_idx, :, 1);    
        err_high = prc_log_r(c_idx, :, 2) - y_val;    
        errorbar(1:cfg.m, y_val, err_low, err_high, '-o', 'Color', colors_A{i}, 'LineWidth', 2, ...
            'MarkerFaceColor', colors_A{i}, 'DisplayName', clean_name(cfg.group_A{i}));
    end
    grid on; xlabel('Spatial Direction Index', 'FontSize', 16); ylabel('log(r_i)', 'FontSize', 16);
    title(sprintf('Group A (Stimulus Delivered) [%s]', band_name), 'FontSize', 16); legend('Location', 'best'); ylim(y_limits);
    
    subplot(1, 2, 2); hold on;
    yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ceiling (BLT B)');
    colors_B = {[0.46 0.67 0.18], [0.30 0.74 0.93]};
    for i = 1:length(cfg.group_B)
        c_idx = find(strcmp(cfg.all_conds, cfg.group_B{i}));
        if isempty(c_idx), continue; end
        y_val = mean_log_r(c_idx, :);              
        err_low = y_val - prc_log_r(c_idx, :, 1);    
        err_high = prc_log_r(c_idx, :, 2) - y_val;    
        errorbar(1:cfg.m, y_val, err_low, err_high, '-s', 'Color', colors_B{i}, 'LineWidth', 2, ...
            'MarkerFaceColor', colors_B{i}, 'DisplayName', clean_name(cfg.group_B{i}));
    end
    grid on; xlabel('Spatial Direction Index', 'FontSize', 16); ylabel('log(r_i)', 'FontSize', 16);
    title(sprintf('Group B (Omission / Anticipation) [%s]', band_name), 'FontSize', 16); legend('Location', 'best'); ylim(y_limits);

    % ----------------------------------------------------------------------
    %% 7. Summary Table (BLT Ref)
    % ----------------------------------------------------------------------
    fprintf('\n========================================================================================\n');
    fprintf('SUMMARY TABLE: BLT REFERENCE (%s BAND)\n', upper(band_name));
    fprintf('----------------------------------------------------------------------------------------\n');
    fprintf('%-18s | %-10s | %-10s | %-18s | %-15s\n', 'Condition', 'Gain (G)', 'log(r_1)', '\Delta Latency (\xi_1)', 'Sustained % (\xi_1)');
    fprintf('----------------------------------------------------------------------------------------\n');
    for c = 1:length(cfg.all_conds)
        cond = cfg.all_conds{c};
        med_G   = median(G_splits(c, :)); med_lr1 = median(log_r_splits(c, 1, :));
        med_lat = median(lat_splits(c, :)); med_sus = median(sust_splits(c, 1, :)) * 100;
        fprintf('%-18s | %-10.3f | %-10.3f | %+7.1f ms         | %6.2f%%\n', cond, med_G, med_lr1, med_lat, med_sus);
    end

    % ----------------------------------------------------------------------
    %% 8. CHECKPOINT 1: Multi-Subject Eigenbasis Extraction (P1 Ref)
    % ----------------------------------------------------------------------
    num_subjects_P1 = length(current_data.P1);
    Xi_all_subjs_P1 = nan(cfg.num_ch, cfg.m, num_subjects_P1);
    lam_frac_all_P1 = nan(cfg.m, num_subjects_P1);
    gap_all_P1      = nan(cfg.m - 1, num_subjects_P1);
    
    for s = 1:num_subjects_P1
        raw_P1 = current_data.P1{s};
        n_trials = size(raw_P1, 3);
        if n_trials < (cfg.nTrialMatch * 2)
            error('Subject %d has insufficient P1 trials (%d < 60)', s, n_trials);
        end
        
        rng(cfg.seed + s);
        perm = randperm(n_trials);
        iA = perm(1:cfg.nTrialMatch);
        
        erpA = mean(raw_P1(:, :, iA), 3, 'omitnan');
        bA   = mean(erpA(:, iBase), 2);
        
        C_A_cent = covCentered(erpA(:, iWin), bA, cfg.alpha);
        
        [V, D] = eig(C_A_cent, 'vector');
        [lam, ord] = sort(D, 'descend');
        V = V(:, ord);
        
        Xi_all_subjs_P1(:, :, s) = V(:, 1:cfg.m);
        lam_frac_all_P1(:, s)   = lam(1:cfg.m) / sum(lam);
        gap_all_P1(:, s)        = (lam(1:cfg.m-1) - lam(2:cfg.m)) ./ lam(1:cfg.m-1);
    end
    
    Xi_aligned_P1 = Xi_all_subjs_P1;
    for s = 2:num_subjects_P1
        for k = 1:cfg.m
            if dot(Xi_aligned_P1(:, k, s), Xi_aligned_P1(:, k, 1)) < 0
                Xi_aligned_P1(:, k, s) = -Xi_aligned_P1(:, k, s);
            end
        end
    end
    Xi_bar_P1 = mean(Xi_aligned_P1, 3);
    
    % --- Deliverable A (P1) ---
    figA_P1 = figure('Position', [100, 100, 1000, 400], 'Name', sprintf('Chk 1 (P1) - A (%s)', band_name));
    subplot(1, 2, 1);
    plot(1:cfg.m, lam_frac_all_P1, '-o', 'LineWidth', 1.5);
    grid on; xlabel('Direction Index', 'FontSize', 16); ylabel('Fraction of Total Variance', 'FontSize', 16);
    title(sprintf('Deliverable A (P1 Ref) [%s]', band_name), 'FontSize', 16);
    
    subplot(1, 2, 2);
    bar(mean(gap_all_P1, 2)); hold on;
    errorbar(1:cfg.m-1, mean(gap_all_P1, 2), std(gap_all_P1, 0, 2), 'k.', 'LineWidth', 1.2);
    grid on; xlabel('Gap Index', 'FontSize', 16); ylabel('Relative Gap', 'FontSize', 16);
    title(sprintf('Mean Subspace Gaps (P1) [%s]', band_name), 'FontSize', 16);

    % --- Deliverable B (P1) ---
    n_rows_P1 = min(num_subjects_P1, 4);
    figure('Position', [50, 50, fig_w, fig_h], 'Name', sprintf('Chk 1 (P1) - B (%s)', band_name));
    for s = 1:n_rows_P1
        for k = 1:cfg.m
            ax = subplot(n_rows_P1, cfg.m, (s-1)*cfg.m + k);
            pos = get(ax, 'Position');
            set(ax, 'Position', [pos(1) * 0.90 + 0.03, pos(2) * 0.95, pos(3) * 0.90, pos(4) * 0.95]);
            topoplot(Xi_aligned_P1(:, k, s), EEG.chanlocs, 'numcontour', 0);
            clim([-0.4 0.4]); colormap('jet');
            if s == 1
                text(0.5, 1.18, sprintf('\\xi_%d', k), 'Units', 'normalized', 'FontSize', 22, ...
                    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            end
            if k == 1
                text(-0.28, 0.5, sprintf('Subj %d', s), 'Units', 'normalized', 'FontSize', 18, ...
                    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
            end
        end
    end
    cb2 = colorbar('Position', [0.88, 0.12, 0.022, 0.72]); 
    cb2.FontSize = 13; cb2.FontWeight = 'bold'; cb2.Ticks = [-0.4, -0.2, 0, 0.2, 0.4];
    cb2.Label.String = 'Spatial Weight (a.u.)'; cb2.Label.FontSize = 15; cb2.Label.FontWeight = 'bold';

    % ----------------------------------------------------------------------
    %% 9. PROJECTIONS & METRICS (Using P1 as Reference)
    % ----------------------------------------------------------------------
    cfg.ref_cond  = 'P1';
    cfg.group_A   = {'P2_500', 'P3_500'};              
    cfg.group_B   = {'P2_2000', 'P3_missing'};         
    cfg.all_conds = [{'P1'}, cfg.group_A, cfg.group_B];

    target_subj = 1; 
    Xi_P1 = Xi_aligned_P1(:, :, target_subj); 

    G_splits_P1       = nan(length(cfg.all_conds), cfg.nRep);
    log_r_splits_P1   = nan(length(cfg.all_conds), cfg.m, cfg.nRep);
    lat_splits_P1     = nan(length(cfg.all_conds), cfg.nRep);
    sust_splits_P1    = nan(length(cfg.all_conds), cfg.m, cfg.nRep);

    for rep = 1:cfg.nRep
        rng(cfg.seed + rep);
        nP1 = size(current_data.P1{target_subj}, 3);
        perm = randperm(nP1);
        iB = perm(cfg.nTrialMatch + (1:cfg.nTrialMatch));
        
        erpB = mean(current_data.P1{target_subj}(:, :, iB), 3, 'omitnan');
        bB   = mean(erpB(:, iBase), 2);
        XB   = erpB(:, iWin) - bB;
        
        v_ceil = sum((Xi_P1' * XB).^2, 2) / T_win;
        V_ceil = sum(v_ceil);
        
        proj_B = Xi_P1(:, 1)' * XB;
        [~, max_idx_B] = max(abs(proj_B));
        t_lat_ceil = t_eval(max_idx_B);
        
        for c = 1:length(cfg.all_conds)
            cond = cfg.all_conds{c};
            if ~isfield(current_data, cond), continue; end
            
            raw_cond = current_data.(cond){target_subj};
            if strcmp(cond, 'P1')
                Xk = XB; bk = bB;
            else
                idx_k = randperm(size(raw_cond, 3), cfg.nTrialMatch);
                erp_k = mean(raw_cond(:, :, idx_k), 3, 'omitnan');
                bk    = mean(erp_k(:, iBase), 2);
                Xk    = erp_k(:, iWin) - bk;
            end
            
            v_k = sum((Xi_P1' * Xk).^2, 2) / T_win;
            V_k = sum(v_k);
            
            G_splits_P1(c, rep)        = V_k / V_ceil;
            log_r_splits_P1(c, :, rep) = log(v_k ./ v_ceil);
            
            proj_k = Xi_P1(:, 1)' * Xk;
            [~, max_idx_k] = max(abs(proj_k));
            lat_splits_P1(c, rep) = (t_eval(max_idx_k) - t_lat_ceil) * 1000; 
            
            xbar = mean(Xk, 2); Xd = Xk - xbar;
            v_sust = (Xi_P1' * xbar).^2;
            v_dyn  = sum((Xi_P1' * Xd).^2, 2) / T_win;
            sust_splits_P1(c, :, rep) = v_sust ./ (v_sust + v_dyn + eps);
        end
    end

    % ----------------------------------------------------------------------
    %% 10. CHECKPOINT 2: Results Display (P1 Ref)
    % ----------------------------------------------------------------------
    mean_log_r_P1 = mean(log_r_splits_P1, 3, 'omitnan'); 
    prc_log_r_P1  = prctile(log_r_splits_P1, [5 95], 3); 
    
    all_low_P1  = prc_log_r_P1(:, :, 1);
    all_high_P1 = prc_log_r_P1(:, :, 2);
    y_min_P1    = min(-1.5, floor(min(all_low_P1(:)) * 1.15));
    y_max_P1    = max( 1.5, ceil(max(all_high_P1(:)) * 1.15));
    y_limits_P1 = [y_min_P1, y_max_P1]; 
    
    figure('Position', [100, 100, 1200, 500], 'Name', sprintf('Chk 2 (P1 Ref) - %s', band_name));
    
    subplot(1, 2, 1); hold on;
    yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ceiling (P1 B)');
    colors_A = {[0.92 0.69 0.12], [0.49 0.18 0.55]}; 
    
    for i = 1:length(cfg.group_A)
        c_idx = find(strcmp(cfg.all_conds, cfg.group_A{i}));
        if isempty(c_idx), continue; end
        y_val = mean_log_r_P1(c_idx, :);              
        err_low = y_val - prc_log_r_P1(c_idx, :, 1);    
        err_high = prc_log_r_P1(c_idx, :, 2) - y_val;    
        errorbar(1:cfg.m, y_val, err_low, err_high, '-o', 'Color', colors_A{i}, 'LineWidth', 2, ...
            'MarkerFaceColor', colors_A{i}, 'DisplayName', clean_name(cfg.group_A{i}));
    end
    grid on; xlabel('Spatial Direction Index', 'FontSize', 16); ylabel('log(r_i)', 'FontSize', 16);
    title(sprintf('Group A vs P1 (Stim Delivered) [%s]', band_name), 'FontSize', 16); legend('Location', 'best'); ylim(y_limits_P1);
    
    subplot(1, 2, 2); hold on;
    yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ceiling (P1 B)');
    colors_B = {[0.46 0.67 0.18], [0.30 0.74 0.93]};
    
    for i = 1:length(cfg.group_B)
        c_idx = find(strcmp(cfg.all_conds, cfg.group_B{i}));
        if isempty(c_idx), continue; end
        y_val = mean_log_r_P1(c_idx, :);              
        err_low = y_val - prc_log_r_P1(c_idx, :, 1);    
        err_high = prc_log_r_P1(c_idx, :, 2) - y_val;    
        errorbar(1:cfg.m, y_val, err_low, err_high, '-s', 'Color', colors_B{i}, 'LineWidth', 2, ...
            'MarkerFaceColor', colors_B{i}, 'DisplayName', clean_name(cfg.group_B{i}));
    end
    grid on; xlabel('Spatial Direction Index', 'FontSize', 16); ylabel('log(r_i)', 'FontSize', 16);
    title(sprintf('Group B vs P1 (Omission) [%s]', band_name), 'FontSize', 16); legend('Location', 'best'); ylim(y_limits_P1);

    % ----------------------------------------------------------------------
    %% 11. Summary Table (P1 Ref)
    % ----------------------------------------------------------------------
    fprintf('\n========================================================================================\n');
    fprintf('SUMMARY TABLE: P1 REFERENCE (%s BAND)\n', upper(band_name));
    fprintf('----------------------------------------------------------------------------------------\n');
    fprintf('%-18s | %-10s | %-10s | %-18s | %-15s\n', 'Condition', 'Gain (G)', 'log(r_1)', '\Delta Latency (\xi_1)', 'Sustained % (\xi_1)');
    fprintf('----------------------------------------------------------------------------------------\n');
    for c = 1:length(cfg.all_conds)
        cond = cfg.all_conds{c};
        med_G   = median(G_splits_P1(c, :)); med_lr1 = median(log_r_splits_P1(c, 1, :));
        med_lat = median(lat_splits_P1(c, :)); med_sus = median(sust_splits_P1(c, 1, :)) * 100;
        fprintf('%-18s | %-10.3f | %-10.3f | %+7.1f ms         | %6.2f%%\n', cond, med_G, med_lr1, med_lat, med_sus);
    end

end % END MASTER FREQUENCY BAND LOOP

disp('All analyses across Raw, Alpha, and Beta bands completed successfully.');