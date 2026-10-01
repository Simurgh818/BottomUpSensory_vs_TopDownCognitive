%% =========================================================================
%% INVESTIGATE COHORT B OUTLIER: MULTI-START K-MEANS & ARTIFACT DIAGNOSTICS
%% =========================================================================
clear; clc; close all;

% 1. Load Precomputed Results from Raw Pipeline
if exist('I:\', 'dir')
    base_output_path = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper\';
elseif exist('H:\', 'dir')
    base_output_path = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper\';
elseif exist('G:\', 'dir')
    base_output_path = 'C:\Users\sdabiri\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper\';
else
    base_output_path = './';
end

raw_results_file = fullfile(base_output_path, 'ch-ch_cov', 'group', 'Raw', 'P1', 'results.mat');
if ~exist(raw_results_file, 'file')
    error('Could not find results.mat at: %s. Run condition_subspace_cov.m first.', raw_results_file);
end

load(raw_results_file, 'G_splits'); % [nCond x nRep x num_subjects]

conditions = {'P1', 'P2_500', 'P3_500', 'P2_2000', 'P3_missing'};
cond_labels = {'Cued', 'Unpred 500', 'Rand 500', 'Unpred 2000', 'Rand Null'};
idx_p3_500 = find(strcmp(conditions, 'P3_500'));
idx_p3_mis = find(strcmp(conditions, 'P3_missing'));
idx_p1     = find(strcmp(conditions, 'P1'));

[nCond, nRep, num_subjs] = size(G_splits);

% Subject Mean Gain across Repetitions: [num_subjs x nCond]
subj_G_matrix = squeeze(mean(G_splits, 2, 'omitnan'))';

% Define Established Cohorts
cohort_A_idx = intersect(4:11, 1:num_subjs);
cohort_B_idx = setdiff(1:num_subjs, cohort_A_idx);

fprintf('Loaded Raw Gain Data for %d subjects across %d repetitions.\n', num_subjs, nRep);
fprintf('Cohort A (Subjs 4-11): n = %d\n', length(cohort_A_idx));
fprintf('Cohort B (Other Subjs): n = %d\n\n', length(cohort_B_idx));

%% =========================================================================
%% 1. PRINCIPLED K-MEANS: MULTI-SEED STABILITY & k=2 vs. k=3
%% =========================================================================
% Feature Space: Gain across the two critical uncertainty conditions
% Feature 1: Rand 500 Gain | Feature 2: Rand Null Gain
X_features = [subj_G_matrix(:, idx_p3_500), subj_G_matrix(:, idx_p3_mis)];

n_seeds = 100;
k_candidates = [2, 3];
seeds = 1000 + (1:n_seeds);

coassign_matrix_k2 = zeros(num_subjs, num_subjs);
coassign_matrix_k3 = zeros(num_subjs, num_subjs);

cluster_assignments_k2 = zeros(num_subjs, n_seeds);
cluster_assignments_k3 = zeros(num_subjs, n_seeds);

sil_scores_k2 = zeros(n_seeds, 1);
sil_scores_k3 = zeros(n_seeds, 1);
inertia_k2    = zeros(n_seeds, 1);
inertia_k3    = zeros(n_seeds, 1);

for s_idx = 1:n_seeds
    rng(seeds(s_idx));
    
    % k = 2
    [idx2, ~, sumd2] = kmeans(X_features, 2, 'Replicates', 10, 'Distance', 'sqeuclidean');
    cluster_assignments_k2(:, s_idx) = idx2;
    sil2 = silhouette(X_features, idx2, 'sqeuclidean');
    sil_scores_k2(s_idx) = mean(sil2);
    inertia_k2(s_idx) = sum(sumd2);
    coassign_matrix_k2 = coassign_matrix_k2 + (idx2 == idx2');
    
    % k = 3
    [idx3, ~, sumd3] = kmeans(X_features, 3, 'Replicates', 10, 'Distance', 'sqeuclidean');
    cluster_assignments_k3(:, s_idx) = idx3;
    sil3 = silhouette(X_features, idx3, 'sqeuclidean');
    sil_scores_k3(s_idx) = mean(sil3);
    inertia_k3(s_idx) = sum(sumd3);
    coassign_matrix_k3 = coassign_matrix_k3 + (idx3 == idx3');
end

coassign_matrix_k2 = coassign_matrix_k2 / n_seeds;
coassign_matrix_k3 = coassign_matrix_k3 / n_seeds;

% Final Representative Clustering (using first seed)
rng(seeds(1));
[final_k2, C_k2] = kmeans(X_features, 2, 'Replicates', 50);
[final_k3, C_k3] = kmeans(X_features, 3, 'Replicates', 50);

fprintf('========================================================================\n');
fprintf('K-MEANS CLUSTERING EVALUATION (Average across %d Random Seeds)\n', n_seeds);
fprintf('========================================================================\n');
fprintf('Model    | Mean Silhouette | Mean Inertia (Sum Sq Dist) | Stability (Co-assign Var)\n');
fprintf('------------------------------------------------------------------------\n');
fprintf('k = 2    | %15.4f | %26.2f | %22.4f\n', mean(sil_scores_k2), mean(inertia_k2), var(coassign_matrix_k2(:)));
fprintf('k = 3    | %15.4f | %26.2f | %22.4f\n', mean(sil_scores_k3), mean(inertia_k3), var(coassign_matrix_k3(:)));
fprintf('------------------------------------------------------------------------\n\n');

% Print Membership Tables
fprintf('CLUSTER MEMBERSHIP BREAKDOWN:\n');
for k_val = 1:2
    fprintf('  k=2, Cluster %d: Subjects %s\n', k_val, mat2str(find(final_k2 == k_val)'));
end
for k_val = 1:3
    fprintf('  k=3, Cluster %d: Subjects %s\n', k_val, mat2str(find(final_k3 == k_val)'));
end
fprintf('\n');

%% =========================================================================
%% 2. OUTLIER IDENTIFICATION & CAUSAL ROOT INSPECTION
%% =========================================================================
% Detect top outlier participant based on distance from median in uncertainty conditions
dist_from_median = sqrt((X_features(:,1) - median(X_features(:,1))).^2 + ...
                        (X_features(:,2) - median(X_features(:,2))).^2);
[max_dist, outlier_subj] = max(dist_from_median);

fprintf('========================================================================\n');
fprintf('OUTLIER INSPECTION: SUBJECT %d (Cohort %s)\n', outlier_subj, ...
    char(categorical(ismember(outlier_subj, cohort_A_idx), [true, false], {'A', 'B'})));
fprintf('========================================================================\n');
fprintf('Raw Gain (Rand 500)  : %8.3f (Cohort B Mean: %.3f | Cohort A Mean: %.3f)\n', ...
    subj_G_matrix(outlier_subj, idx_p3_500), mean(subj_G_matrix(cohort_B_idx, idx_p3_500)), mean(subj_G_matrix(cohort_A_idx, idx_p3_500)));
fprintf('Raw Gain (Rand Null) : %8.3f (Cohort B Mean: %.3f | Cohort A Mean: %.3f)\n', ...
    subj_G_matrix(outlier_subj, idx_p3_mis), mean(subj_G_matrix(cohort_B_idx, idx_p3_mis)), mean(subj_G_matrix(cohort_A_idx, idx_p3_mis)));

% Recalculate Cohort B Mean WITHOUT the outlier
cohort_B_ex = setdiff(cohort_B_idx, outlier_subj);
fprintf('\nIMPACT OF REMOVING OUTLIER (Subj %d) FROM COHORT B:\n', outlier_subj);
fprintf('  Rand 500 Gain  : %.3f  -->  %.3f (Drop of %.1f%%)\n', ...
    mean(subj_G_matrix(cohort_B_idx, idx_p3_500)), mean(subj_G_matrix(cohort_B_ex, idx_p3_500)), ...
    100*(mean(subj_G_matrix(cohort_B_idx, idx_p3_500)) - mean(subj_G_matrix(cohort_B_ex, idx_p3_500))) / mean(subj_G_matrix(cohort_B_idx, idx_p3_500)));
fprintf('  Rand Null Gain : %.3f  -->  %.3f (Drop of %.1f%%)\n\n', ...
    mean(subj_G_matrix(cohort_B_idx, idx_p3_mis)), mean(subj_G_matrix(cohort_B_ex, idx_p3_mis)), ...
    100*(mean(subj_G_matrix(cohort_B_idx, idx_p3_mis)) - mean(subj_G_matrix(cohort_B_ex, idx_p3_mis))) / mean(subj_G_matrix(cohort_B_idx, idx_p3_mis)));

%% =========================================================================
%% 3. FIGURE 1: K-MEANS CLUSTERING & MODEL SELECTION
%% =========================================================================
fig_km = figure('Position', [100, 100, 1600, 500], 'Name', 'K-Means Comparison');
t_km = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'normal');

% Tile 1: 2D Feature Space with k=2 Solution
nexttile; hold on; set(gca, 'FontSize', 16);
colors_k2 = {[0.2 0.6 0.2], [0.8 0.2 0.2]};
for k_val = 1:2
    members = find(final_k2 == k_val);
    scatter(X_features(members, 1), X_features(members, 2), 120, colors_k2{k_val}, 'filled', ...
        'DisplayName', sprintf('Cluster %d (n=%d)', k_val, length(members)));
end
for s = 1:num_subjs
    text(X_features(s, 1) + 0.3, X_features(s, 2), sprintf('S%d', s), 'FontSize', 12, 'FontWeight', 'bold');
end
xlabel('Raw Gain: Rand 500 (P3)'); ylabel('Raw Gain: Rand Null (P3)');
title('k = 2 Partition', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'northwest'); grid on;

% Tile 2: 2D Feature Space with k=3 Solution
nexttile; hold on; set(gca, 'FontSize', 16);
colors_k3 = {[0.2 0.6 0.2], [0.1 0.4 0.8], [0.8 0.2 0.2]};
for k_val = 1:3
    members = find(final_k3 == k_val);
    scatter(X_features(members, 1), X_features(members, 2), 120, colors_k3{k_val}, 'filled', ...
        'DisplayName', sprintf('Cluster %d (n=%d)', k_val, length(members)));
end
for s = 1:num_subjs
    text(X_features(s, 1) + 0.3, X_features(s, 2), sprintf('S%d', s), 'FontSize', 12, 'FontWeight', 'bold');
end
xlabel('Raw Gain: Rand 500 (P3)'); ylabel('Raw Gain: Rand Null (P3)');
title('k = 3 Partition (Singleton Isolation)', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'northwest'); grid on;

% Tile 3: Co-assignment Matrix (Stability across 100 seeds for k=3)
nexttile;
imagesc(coassign_matrix_k3); colormap('hot'); colorbar; clim([0 1]);
set(gca, 'FontSize', 16, 'XTick', 1:num_subjs, 'YTick', 1:num_subjs);
xlabel('Subject ID'); ylabel('Subject ID');
title('k = 3 Co-Assignment Matrix (100 Seeds)', 'FontSize', 18, 'FontWeight', 'bold');
axis square;

out_km_file = fullfile(base_output_path, 'ch-ch_cov', 'Figure_KMeans_Participant_Clustering');
saveas(fig_km, [out_km_file, '.png']); saveas(fig_km, [out_km_file, '.svg']);
fprintf('Saved k-Means evaluation to: %s\n', out_km_file);

%% =========================================================================
%% 4. FIGURE 2: PARTICIPANT-BY-PARTICIPANT GAIN PROFILES
%% =========================================================================
fig_diag = figure('Position', [100, 100, 1600, 600], 'Name', 'Participant Gain Profiles');
t_diag = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'normal');

% Panel A: Rand 500 Gain per Participant
nexttile; hold on; set(gca, 'FontSize', 16);
yline(1.0, 'k--', 'LineWidth', 2, 'DisplayName', 'Cued Ceiling');
bar_cols = repmat([0.2 0.4 0.8], num_subjs, 1);
bar_cols(cohort_A_idx, :) = repmat([0.5 0.5 0.5], length(cohort_A_idx), 1);
bar_cols(outlier_subj, :) = [0.85 0.2 0.2]; % Highlight outlier

for s = 1:num_subjs
    b = bar(s, subj_G_matrix(s, idx_p3_500), 0.7);
    b.FaceColor = bar_cols(s, :);
end
xticks(1:num_subjs);
xlabel('Participant ID'); ylabel('Global Gain (G)');
title('Rand 500 (P3) Gain per Participant', 'FontSize', 18, 'FontWeight', 'bold');
grid on;

% Panel B: Rand Null Gain per Participant
nexttile; hold on; set(gca, 'FontSize', 16);
yline(1.0, 'k--', 'LineWidth', 2, 'DisplayName', 'Cued Ceiling');
for s = 1:num_subjs
    b = bar(s, subj_G_matrix(s, idx_p3_mis), 0.7);
    b.FaceColor = bar_cols(s, :);
end
xticks(1:num_subjs);
xlabel('Participant ID'); ylabel('Global Gain (G)');
title('Rand Null (P3 Missing) Gain per Participant', 'FontSize', 18, 'FontWeight', 'bold');
grid on;

sgtitle('Participant-Level Gain Decomposition: Gray = Cohort A, Blue = Cohort B, Red = Driver', ...
    'FontSize', 20, 'FontWeight', 'bold');

out_diag_file = fullfile(base_output_path, 'ch-ch_cov', 'Figure_Participant_Gain_Decomposition');
saveas(fig_diag, [out_diag_file, '.png']); saveas(fig_diag, [out_diag_file, '.svg']);
fprintf('Saved Participant Gain Decomposition to: %s\n', out_diag_file);