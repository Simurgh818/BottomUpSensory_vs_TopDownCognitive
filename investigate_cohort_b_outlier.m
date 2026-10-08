%% =========================================================================
%% OUTLIER HUNTER: K-MEANS DISCOVERY & 3-HYPOTHESIS DIAGNOSTIC
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
cohort_B_idx = setdiff(1:num_subjs, [cohort_A_idx,14]);

%% =========================================================================
%% 3. IDENTIFY THE NEW OUTLIER (IGNORING EXCLUDED SUBJECTS)
%% =========================================================================
% Define the pool of valid active subjects
valid_subjs = [cohort_A_idx, cohort_B_idx];
excluded_subjs = setdiff(1:num_subjs, valid_subjs);

X_features = [subj_G_matrix(:, idx_p3_500), subj_G_matrix(:, idx_p3_mis)];

% Calculate distance from the median of ONLY the valid subjects
dist_from_median = sqrt((X_features(:,1) - median(X_features(valid_subjs, 1))).^2 + ...
                        (X_features(:,2) - median(X_features(valid_subjs, 2))).^2);

% Force excluded subjects (e.g., S14) to have a distance of -1 so they are never picked
dist_from_median(excluded_subjs) = -1;

[sorted_dist, sort_idx] = sort(dist_from_median, 'descend');
outlier_subj = sort_idx(1); % Automatically grab the absolute worst valid offender

fprintf('========================================================================\n');
fprintf('OUTLIER DISCOVERY REPORT\n');
fprintf('========================================================================\n');
fprintf('Top 3 Outliers driving variance (Excluding Subject 14):\n');
for i = 1:3
    s_id = sort_idx(i);
    grp = 'B'; if ismember(s_id, cohort_A_idx), grp = 'A'; end
    fprintf('  #%d: Subject %02d (Cohort %s) | Rand 500 Gain: %6.2f | Rand Null Gain: %6.2f\n', ...
        i, s_id, grp, subj_G_matrix(s_id, idx_p3_500), subj_G_matrix(s_id, idx_p3_mis));
end
fprintf('\n=> Running deep diagnostic on worst active offender: Subject %d\n', outlier_subj);

%% =========================================================================
%% 4. PLOT PARTICIPANT GAIN PROFILES (BAR CHARTS - FIXED Y-AXIS)
%% =========================================================================
fig_diag = figure('Position', [100, 100, 1600, 600], 'Name', 'Participant Gain Profiles');
t_diag = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'normal');

bar_cols = repmat([0.2 0.4 0.8], num_subjs, 1); % Default Cohort B (Blue)
bar_cols(cohort_A_idx, :) = repmat([0.5 0.5 0.5], length(cohort_A_idx), 1); % Cohort A (Gray)

% Mask excluded subjects (e.g., S14) so they don't break the Y-axis scale
plot_G_matrix = subj_G_matrix;
if ~isempty(excluded_subjs)
    plot_G_matrix(excluded_subjs, :) = 0; 
    bar_cols(excluded_subjs, :) = repmat([1 1 1], length(excluded_subjs), 1); % Invisible
end
bar_cols(outlier_subj, :) = [0.85 0.2 0.2]; % Highlight top active outlier (Red)

% Panel A: Rand 500
nexttile; hold on; set(gca, 'FontSize', 16);
yline(1.0, 'k--', 'LineWidth', 2);
for s = 1:num_subjs, b = bar(s, plot_G_matrix(s, idx_p3_500), 0.7); b.FaceColor = bar_cols(s, :); end
xticks(1:num_subjs); xlabel('Participant ID'); ylabel('Global Gain (G)');
title('Rand 500 (P3) Gain per Participant', 'FontSize', 18, 'FontWeight', 'bold'); grid on;

% Panel B: Rand Null
nexttile; hold on; set(gca, 'FontSize', 16);
yline(1.0, 'k--', 'LineWidth', 2);
for s = 1:num_subjs, b = bar(s, plot_G_matrix(s, idx_p3_mis), 0.7); b.FaceColor = bar_cols(s, :); end
xticks(1:num_subjs); xlabel('Participant ID'); ylabel('Global Gain (G)');
title('Rand Null (P3 Missing) Gain per Participant', 'FontSize', 18, 'FontWeight', 'bold'); grid on;

sgtitle(sprintf('Gain Decomposition: Gray = Cohort A, Blue = Cohort B, Red = Driver (S%d)', outlier_subj), ...
    'FontSize', 20, 'FontWeight', 'bold');

%% =========================================================================
%% 5. LOAD RAW DATA FOR DIAGNOSTIC (ALL CONDITIONS FOR LINE PLOT)
%% =========================================================================
fprintf('\nLoading raw EEGLAB data for diagnostic...\n');
data_all_conds = struct();
data_all_conds.Raw = struct();

% Added 'P2' to load all conditions for the 5-point line plot
conds_to_load = {'P1', 'P2', 'P3'}; 
for c_idx = 1:length(conds_to_load)
    cond_name = conds_to_load{c_idx};
    in_dir = fullfile(input_path, cond_name);
    set_files = dir(fullfile(in_dir, '*.set'));
    names_sorted = sort(cellstr({set_files.name})');
    
    for s = 1:num_subjs
        file_to_load = names_sorted{s};
        clean_file_check = fullfile(input_path, [cond_name, '_cleaned'], file_to_load);
        
        % Check if an ICA-cleaned file exists
        if exist(clean_file_check, 'file')
            fprintf('Loading ICA-CLEANED Subj %d/%d (%s) for condition: %s\n', ...
                s, num_subjs, file_to_load, cond_name);
            EEG = pop_loadset('filename', file_to_load, 'filepath', fullfile(input_path, [cond_name, '_cleaned']));
        else
            fprintf('Loading Subj %d/%d (%s) for condition: %s\n', ...
                s, num_subjs, file_to_load, cond_name);
            EEG = pop_loadset('filename', file_to_load, 'filepath', in_dir);
        end
        
        if strcmp(cond_name, 'P1')
            data_all_conds.Raw.P1{s} = EEG.data;
        elseif strcmp(cond_name, 'P2')
            data_all_conds.Raw.P2_500{s}     = EEG.data(:, :, 1:2:EEG.trials);
            data_all_conds.Raw.P2_2000{s}    = EEG.data(:, :, 2:2:EEG.trials);
        elseif strcmp(cond_name, 'P3')
            data_all_conds.Raw.P3_500{s}     = EEG.data(:, :, 1:2:EEG.trials);
            data_all_conds.Raw.P3_missing{s} = EEG.data(:, :, 2:2:EEG.trials);
        end
    end
end

% Robust time reconstruction
nPnts = size(data_all_conds.Raw.P1{1}, 2);
time_s = -1.000 + (0:(nPnts - 1)) / 500; 

%% =========================================================================
%% 6. CALCULATE DIAGNOSTIC METRICS (TARGETING P3_500 & ALL BASELINES)
%% =========================================================================
cfg = struct('win', [0 1.0], 'base', [-1.0 -0.8], 'm', 6);
iWin  = time_s >= cfg.win(1)  & time_s < cfg.win(2);
iBase = time_s >= cfg.base(1) & time_s < cfg.base(2);
covCentered = @(X, b) (( (X - b) - mean(X - b, 2) ) * ( (X - b) - mean(X - b, 2) )') / size(X, 2);

V_ceil_all   = zeros(num_subjs, 1);
Base_Var_all = zeros(num_subjs, 5);
cond_keys    = {'P1', 'P2_500', 'P3_500', 'P2_2000', 'P3_missing'};

for s = valid_subjs % ONLY PROCESS VALID SUBJECTS
    raw_ref = data_all_conds.Raw.P1{s};
    half_n = floor(size(raw_ref, 3) / 2);
    erpA = mean(raw_ref(:, :, 1:half_n), 3, 'omitnan');
    erpB = mean(raw_ref(:, :, (half_n+1):(2*half_n)), 3, 'omitnan');
    
    bA = mean(erpA(:, iBase), 2); bB = mean(erpB(:, iBase), 2);
    C_A = covCentered(erpA(:, iWin), bA);
    [V, D] = eig(C_A, 'vector'); [~, ord] = sort(D, 'descend');
    Xi = V(:, ord(1:cfg.m));
    
    XB = erpB(:, iWin) - bB;
    V_ceil_all(s) = sum(sum((Xi' * XB).^2, 2) / nnz(iWin));
    
    % Calculate baseline noise for all 5 conditions for the line plot
    for c = 1:5
        erp_all = mean(data_all_conds.Raw.(cond_keys{c}){s}, 3, 'omitnan');
        Base_Var_all(s, c) = mean(var(erp_all(:, iBase), 0, 2));
    end
end

% FIX: Only use valid subjects to calculate normative medians!
normative_subjs = setdiff(valid_subjs, outlier_subj);
median_Vceil = median(V_ceil_all(normative_subjs));
Vceil_ratio = V_ceil_all(outlier_subj) / median_Vceil;

% Evaluate drift specifically for P3_500 vs P1
base_drift  = Base_Var_all(outlier_subj, 3) / (Base_Var_all(outlier_subj, 1) + eps); 

fprintf('\n========================================================================\n');
fprintf(' DIAGNOSTIC RESULT FOR SUBJECT %d (RAND 500)\n', outlier_subj);
fprintf('========================================================================\n');
fprintf('1. Denominator (V_ceil) Ratio : %6.2f%% of Normal\n', Vceil_ratio * 100);
fprintf('2. Baseline Noise Drift (P3 500 / P1): %6.2f x\n', base_drift);

if Vceil_ratio < 0.25 && base_drift <= 3.0
    fprintf('>>> VERDICT: DENOMINATOR COLLAPSE (Artifactual Math)\n');
elseif base_drift > 3.0
    fprintf('>>> VERDICT: NON-NEURAL VOLTAGE ARTIFACT (Noise Drift)\n');
else
    fprintf('>>> VERDICT: GENUINE HYPER-RESPONDER\n');
end

%% =========================================================================
%% 7. PLOT 3-PANEL DIAGNOSTIC (LINE PLOT)
%% =========================================================================
fig_3panel = figure('Position', [100, 100, 1600, 500], 'Name', sprintf('S%d Diagnostic', outlier_subj));
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'normal');

% Panel A: V_ceil
nexttile; hold on; set(gca, 'FontSize', 14);
bar(normative_subjs, V_ceil_all(normative_subjs), 'FaceColor', [0.6 0.6 0.6]);
bar(outlier_subj, V_ceil_all(outlier_subj), 'FaceColor', [0.1 0.4 0.8]); % Blue for Target
yline(median_Vceil, 'k--', 'LineWidth', 2);
xlabel('Participant ID'); ylabel('V_{ceil} (Projected \muV^2)');
title('Reference Denominator (V_{ceil})', 'FontSize', 16, 'FontWeight', 'bold'); grid on;
xlim([0, num_subjs+1]);

% Panel B: Noise Drift (Line Plot)
nexttile; hold on; set(gca, 'FontSize', 14);
plot(1:5, mean(Base_Var_all(normative_subjs, :), 1), 'k-o', 'LineWidth', 2, 'DisplayName', 'Cohort Mean');
plot(1:5, Base_Var_all(outlier_subj, :), 'b-s', 'LineWidth', 2.5, 'MarkerSize', 8, 'DisplayName', sprintf('Subject %d', outlier_subj));
xticks(1:5); xticklabels({'P1', 'P2_{500}', 'P3_{500}', 'P2_{2000}', 'P3_{mis}'});
ylabel('Pre-Cue Baseline Variance (\muV^2)');
title('Pre-Cue Noise Stability', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best'); grid on;

% Panel C: Waveform Check
chan_idx = 13; % Default to C3/Cz equivalent
erp_P1     = squeeze(mean(data_all_conds.Raw.P1{outlier_subj}(chan_idx, :, :), 3, 'omitnan'));
erp_P3_500 = squeeze(mean(data_all_conds.Raw.P3_500{outlier_subj}(chan_idx, :, :), 3, 'omitnan'));

nexttile; hold on; set(gca, 'FontSize', 14); xline(0, 'k--');
plot(time_s(:), erp_P1(:), 'k-', 'LineWidth', 2, 'DisplayName', 'Cued (P1)');
plot(time_s(:), erp_P3_500(:), 'Color', [0.9 0.6 0.1], 'LineWidth', 2, 'DisplayName', 'Rand 500'); 
xlim([-0.2, 0.8]); xlabel('Time (s)'); ylabel('Voltage (\muV)');
title(sprintf('Subject %d Waveforms (Rand 500)', outlier_subj), 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best'); grid on;

sgtitle(sprintf('Subject %d Diagnostic Profile', outlier_subj), 'FontSize', 18, 'FontWeight', 'bold');