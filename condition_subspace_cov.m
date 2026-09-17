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

%% 0.2 Apply Alpha and Beta Zero-Phase Filtering
disp('Applying Zero-Phase FIR Filters for Alpha (8-12 Hz) and Beta (13-30 Hz)...');

fn = fs / 2; 
ord_alpha = round(0.250 * fs); % 250 ms
b_alpha   = fir1(ord_alpha, [8 12] / fn, 'bandpass');

ord_beta  = round(0.125 * fs); % 125 ms
b_beta    = fir1(ord_beta, [13 30] / fn, 'bandpass');

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
        data_all_conds.Alpha.(c_name){s} = permute(reshape(X_alpha_2D, nt, nc, ntr), [2, 1, 3]);
        
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
%% MASTER LOOP: PROCESS FREQUENCY BANDS AND BASELINES (BLT & P1)
%% =========================================================================
bands_to_process = {'Raw', 'Alpha', 'Beta'};
refs_to_process  = {'BLT', 'P1'};

for band_idx = 1:length(bands_to_process)
    band_name = bands_to_process{band_idx};
    current_data = data_all_conds.(band_name); 
    num_subjects = length(current_data.BLT);
    
    for ref_idx = 1:length(refs_to_process)
        cfg.ref_cond = refs_to_process{ref_idx};
        
        fprintf('\n\n======================================================\n');
        fprintf('   BAND: %s | REFERENCE: %s \n', upper(band_name), cfg.ref_cond);
        fprintf('======================================================\n');
        
        % Set dynamic groups based on reference condition
        if strcmp(cfg.ref_cond, 'BLT')
            cfg.group_A = {'P1', 'P2_500', 'P3_500'};
        else
            cfg.group_A = {'P2_500', 'P3_500'};
        end
        cfg.group_B   = {'P2_2000', 'P3_missing'};
        cfg.all_conds = [{cfg.ref_cond}, cfg.group_A, cfg.group_B];
        
        % Setup Directory Structure[cite: 1]
        group_dir = fullfile(output_path, 'group', band_name, cfg.ref_cond);
        if ~exist(group_dir, 'dir'), mkdir(group_dir); end

        % ------------------------------------------------------------------
        %% CHECKPOINT 1: Group Eigenbasis Evaluation[cite: 1]
        % ------------------------------------------------------------------
        Xi_all_subjs = nan(cfg.num_ch, cfg.m, num_subjects);
        lam_frac_all = nan(cfg.m, num_subjects);
        gap_all      = nan(cfg.m - 1, num_subjects);
        
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
            C_A_cent = covCentered(erpA(:, iWin), bA, cfg.alpha);
            
            [V, D] = eig(C_A_cent, 'vector');
            [lam, ord] = sort(D, 'descend');
            V = V(:, ord);
            
            Xi_all_subjs(:, :, s) = V(:, 1:cfg.m);
            lam_frac_all(:, s)    = lam(1:cfg.m) / sum(lam);
            gap_all(:, s)         = (lam(1:cfg.m-1) - lam(2:cfg.m)) ./ lam(1:cfg.m-1);
        end
        
        % Sign-align eigenvectors to Subject 1
        Xi_aligned = Xi_all_subjs;
        for s = 2:num_subjects
            for k = 1:cfg.m
                if dot(Xi_aligned(:, k, s), Xi_aligned(:, k, 1)) < 0
                    Xi_aligned(:, k, s) = -Xi_aligned(:, k, s);
                end
            end
        end
        Xi_bar = mean(Xi_aligned, 3);
        
        % --- Save Group Level Basis Data ---
        save(fullfile(group_dir, 'spectra.mat'), 'lam_frac_all', 'gap_all');
        
        sim_matrix = nan(num_subjects, cfg.m);
        for s = 1:num_subjects
            for k = 1:cfg.m
                sim_matrix(s, k) = abs(dot(Xi_aligned(:, k, s), Xi_bar(:, k))) / (norm(Xi_aligned(:, k, s)) * norm(Xi_bar(:, k)));
            end
        end
        save(fullfile(group_dir, 'eigvec_similarity.mat'), 'sim_matrix');
        
        % --- Plot & Save Checkpoint 1 Figures ---
        figA = figure('Position', [100, 100, 1000, 400], 'Visible', 'off');
        subplot(1, 2, 1); plot(1:cfg.m, lam_frac_all, '-o', 'LineWidth', 1.5); grid on;
        xlabel('Direction Index'); ylabel('Fraction of Total Variance'); title(sprintf('Eigenvalue Spectra [%s|%s]', band_name, cfg.ref_cond));
        subplot(1, 2, 2); bar(mean(gap_all, 2)); hold on; errorbar(1:cfg.m-1, mean(gap_all, 2), std(gap_all, 0, 2), 'k.', 'LineWidth', 1.2);
        grid on; xlabel('Gap Index'); ylabel('Relative Gap'); title('Mean Subspace Gaps');
        saveas(figA, fullfile(group_dir, 'Deliverable_A_Spectra.png')); close(figA);
        
        n_rows = min(num_subjects, 4);
        figB = figure('Position', [50, 50, 200 * cfg.m + 140, 190 * n_rows + 40], 'Visible', 'off');
        for s = 1:n_rows
            for k = 1:cfg.m
                ax = subplot(n_rows, cfg.m, (s-1)*cfg.m + k);
                pos = get(ax, 'Position'); set(ax, 'Position', [pos(1) * 0.90 + 0.03, pos(2) * 0.95, pos(3) * 0.90, pos(4) * 0.95]);
                topoplot(Xi_aligned(:, k, s), EEG.chanlocs, 'numcontour', 0); clim([-0.4 0.4]); colormap('jet');
                if s == 1, text(0.5, 1.18, sprintf('\\xi_%d', k), 'Units', 'normalized', 'FontSize', 22, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); end
                if k == 1, text(-0.28, 0.5, sprintf('Subj %d', s), 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90); end
            end
        end
        cb = colorbar('Position', [0.88, 0.12, 0.022, 0.72]); cb.FontSize = 13; cb.FontWeight = 'bold'; cb.Ticks = [-0.4, -0.2, 0, 0.2, 0.4]; cb.Label.String = 'Spatial Weight (a.u.)';
        saveas(figB, fullfile(group_dir, 'Deliverable_B_Topos.png')); close(figB);
        
        figC = figure('Position', [150, 150, 600, 400], 'Visible', 'off');
        bar(mean(sim_matrix, 1)); hold on; errorbar(1:cfg.m, mean(sim_matrix, 1), std(sim_matrix, 0, 1), 'k.', 'LineWidth', 1.5);
        grid on; xlabel('Spatial Direction Index'); ylabel('Cosine Similarity'); title('Consistency Across Subjects'); ylim([0 1]);
        saveas(figC, fullfile(group_dir, 'Deliverable_C_Consistency.png')); close(figC);

        % ------------------------------------------------------------------
        %% PROJECTIONS: Run for ALL Subjects & All Reps[cite: 1]
        % ------------------------------------------------------------------
        % Group arrays: [Conditions x Metric/Dims x Reps x Subjects]
        G_splits     = nan(length(cfg.all_conds), cfg.nRep, num_subjects);
        log_r_splits = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
        lat_splits   = nan(length(cfg.all_conds), cfg.nRep, num_subjects);
        sust_splits  = nan(length(cfg.all_conds), cfg.m, cfg.nRep, num_subjects);
        
        for subj = 1:num_subjects
            subj_dir = fullfile(output_path, sprintf('subj%02d', subj), band_name, cfg.ref_cond);
            if ~exist(subj_dir, 'dir'), mkdir(subj_dir); end
            save(fullfile(subj_dir, 'cfg.mat'), 'cfg');
            
            % Stable alignment basis for this subject
            Xi_stable = Xi_aligned(:, :, subj);
            
            for rep = 1:cfg.nRep
                rng(cfg.seed + rep + subj*1000); % Unique split seed
                
                % 1. Dynamic Split-Half
                raw_ref = current_data.(cfg.ref_cond){subj};
                nRef = size(raw_ref, 3);
                perm = randperm(nRef);
                iA = perm(1:cfg.nTrialMatch);
                iB = perm(cfg.nTrialMatch + (1:cfg.nTrialMatch));
                
                % Half A Basis recalculation per rep
                erpA = mean(raw_ref(:, :, iA), 3, 'omitnan');
                bA   = mean(erpA(:, iBase), 2);
                C_A_cent = covCentered(erpA(:, iWin), bA, cfg.alpha);
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
                
                basis_vars = struct('Xi', Xi_rep, 'lam', lam, 'lam_frac', lam_frac_rep, 'gap', gap_rep, 'alpha', cfg.alpha);
                
                % Half B Ceiling Evaluation
                erpB = mean(raw_ref(:, :, iB), 3, 'omitnan');
                bB   = mean(erpB(:, iBase), 2);
                XB   = erpB(:, iWin) - bB;
                
                v_ceil = sum((Xi_rep' * XB).^2, 2) / T_win;
                V_ceil = sum(v_ceil);
                proj_B = Xi_rep(:, 1)' * XB;
                [~, max_idx_B] = max(abs(proj_B));
                t_lat_ceil = t_eval(max_idx_B);
                
                % 2. Projections
                split_vars = struct('iA', iA, 'iB', iB);
                proj_vars  = struct('v_ceil', v_ceil);
                
                for c = 1:length(cfg.all_conds)
                    cond = cfg.all_conds{c};
                    if ~isfield(current_data, cond), continue; end
                    
                    raw_cond = current_data.(cond){subj};
                    if strcmp(cond, cfg.ref_cond)
                        Xk = XB; bk = bB; idx_k = iB;
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
                    
                    % Store in Group Arrays
                    G_splits(c, rep, subj)        = G_k;
                    log_r_splits(c, :, rep, subj) = log_r_k;
                    
                    proj_k = Xi_rep(:, 1)' * Xk;
                    [~, max_idx_k] = max(abs(proj_k));
                    lat_splits(c, rep, subj) = (t_eval(max_idx_k) - t_lat_ceil) * 1000;
                    
                    xbar = mean(Xk, 2); Xd = Xk - xbar;
                    v_sust = (Xi_rep' * xbar).^2;
                    v_dyn  = sum((Xi_rep' * Xd).^2, 2) / T_win;
                    sust_splits(c, :, rep, subj) = v_sust ./ (v_sust + v_dyn + eps);
                end
                
                % Save files per subject per rep[cite: 1]
                rep_str = sprintf('%02d', rep);
                save(fullfile(subj_dir, ['split_rep' rep_str '.mat']), '-struct', 'split_vars');
                save(fullfile(subj_dir, ['basis_rep' rep_str '.mat']), '-struct', 'basis_vars');
                save(fullfile(subj_dir, ['proj_rep' rep_str '.mat']), '-struct', 'proj_vars');
            end
        end

        % ------------------------------------------------------------------
        %% CHECKPOINT 2: Group Level Stats, Tables, and Figures[cite: 1]
        % ------------------------------------------------------------------
        % Average across reps first (within subject), then compute Group Stats
        subj_mean_log_r = squeeze(mean(log_r_splits, 3, 'omitnan')); % [nCond x m x Subjs]
        group_mean_log_r = mean(subj_mean_log_r, 3, 'omitnan');      % [nCond x m]
        group_se_log_r = std(subj_mean_log_r, 0, 3, 'omitnan') ./ sqrt(num_subjects);
        
        figChk2 = figure('Position', [100, 100, 1200, 500], 'Visible', 'off');
        
        % Dynamic Y-Limits
        y_min = min(-1.5, floor(min(group_mean_log_r(:) - group_se_log_r(:)) * 1.15));
        y_max = max( 1.5, ceil(max(group_mean_log_r(:) + group_se_log_r(:)) * 1.15));
        y_limits = [y_min, y_max];
        
        % Group A Plot
        subplot(1, 2, 1); hold on;
        yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('Ceiling (%s B)', cfg.ref_cond));
        colors_A = {[0.85 0.32 0.09], [0.92 0.69 0.12], [0.49 0.18 0.55]};
        for i = 1:length(cfg.group_A)
            c_idx = find(strcmp(cfg.all_conds, cfg.group_A{i}));
            if isempty(c_idx), continue; end
            y_val = group_mean_log_r(c_idx, :);              
            err = group_se_log_r(c_idx, :);    
            col_idx = mod(i-1, 3) + 1;
            errorbar(1:cfg.m, y_val, err, err, '-o', 'Color', colors_A{col_idx}, 'LineWidth', 2, ...
                'MarkerFaceColor', colors_A{col_idx}, 'DisplayName', clean_name(cfg.group_A{i}));
        end
        grid on; xlabel('Spatial Direction Index', 'FontSize', 16); ylabel('log(r_i) \pm SEM', 'FontSize', 16);
        title(sprintf('Group A: Stimulus Delivered (%s %s)', band_name, cfg.ref_cond), 'FontSize', 16); legend('Location', 'best'); ylim(y_limits);
        
        % Group B Plot
        subplot(1, 2, 2); hold on;
        yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('Ceiling (%s B)', cfg.ref_cond));
        colors_B = {[0.46 0.67 0.18], [0.30 0.74 0.93]};
        for i = 1:length(cfg.group_B)
            c_idx = find(strcmp(cfg.all_conds, cfg.group_B{i}));
            if isempty(c_idx), continue; end
            y_val = group_mean_log_r(c_idx, :);              
            err = group_se_log_r(c_idx, :);    
            errorbar(1:cfg.m, y_val, err, err, '-s', 'Color', colors_B{i}, 'LineWidth', 2, ...
                'MarkerFaceColor', colors_B{i}, 'DisplayName', clean_name(cfg.group_B{i}));
        end
        grid on; xlabel('Spatial Direction Index', 'FontSize', 16); ylabel('log(r_i) \pm SEM', 'FontSize', 16);
        title(sprintf('Group B: Omission (%s %s)', band_name, cfg.ref_cond), 'FontSize', 16); legend('Location', 'best'); ylim(y_limits);

        saveas(figChk2, fullfile(group_dir, 'Checkpoint2_Subspace_Redistribution.png')); close(figChk2);

        % --- Export Group Results and Summary Table ---
        save(fullfile(group_dir, 'results.mat'), 'log_r_splits', 'G_splits', 'lat_splits', 'sust_splits');
        
        fid = fopen(fullfile(group_dir, 'summary_table.txt'), 'w');
        for out = [1, fid]
            fprintf(out, '\n========================================================================================\n');
            fprintf(out, 'SUMMARY TABLE: %s REFERENCE (BAND: %s) [Average over %d subjects]\n', cfg.ref_cond, upper(band_name), num_subjects);
            fprintf(out, '----------------------------------------------------------------------------------------\n');
            fprintf(out, '%-18s | %-10s | %-10s | %-18s | %-15s\n', 'Condition', 'Gain (G)', 'log(r_1)', '\Delta Latency (\xi_1)', 'Sustained % (\xi_1)');
            fprintf(out, '----------------------------------------------------------------------------------------\n');
            for c = 1:length(cfg.all_conds)
                cond = cfg.all_conds{c};
                % Average across subjects and reps
                med_G   = mean(G_splits(c, :), 'all', 'omitnan'); 
                med_lr1 = group_mean_log_r(c, 1);
                med_lat = mean(lat_splits(c, :), 'all', 'omitnan'); 
                med_sus = mean(sust_splits(c, 1, :), 'all', 'omitnan') * 100;
                fprintf(out, '%-18s | %-10.3f | %-10.3f | %+7.1f ms         | %6.2f%%\n', cond, med_G, med_lr1, med_lat, med_sus);
            end
        end
        fclose(fid);
        
    end % END REFERENCE LOOP
end % END BAND LOOP

disp('All analyses across Raw, Alpha, and Beta bands completed successfully.');