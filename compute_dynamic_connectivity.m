function [sub_ch_A_all, sub_ch_B_all, sub_pc_A_all, sub_pc_B_all, window_centers, k_opt, pc_labels] = ...
    compute_dynamic_connectivity(trialsA_raw, trialsB_raw, t_win, time_ms_eeg, fs, window_size_ms, step_size_ms, n_perms)

    nA = size(trialsA_raw, 3);
    nB = size(trialsB_raw, 3);
    num_ch = size(trialsA_raw, 1);
    total_trials = nA + nB;

    % 1. GLOBAL Z-SCORE PER TRIAL
    trialsA = zeros(size(trialsA_raw));
    trialsB = zeros(size(trialsB_raw));
    for tr = 1:nA, trialsA(:,:,tr) = zscore(trialsA_raw(:,:,tr), 0, 2); end
    for tr = 1:nB, trialsB(:,:,tr) = zscore(trialsB_raw(:,:,tr), 0, 2); end

    avgA = mean(trialsA, 3, 'omitnan');
    avgB = mean(trialsB, 3, 'omitnan');

    % 2. DIMENSIONALITY REDUCTION (SHARED MANIFOLD via dPCA)
    data_combined = [avgA, avgB];
    data_combined_centered = data_combined - mean(data_combined, 2);
    
    [~, ~, V] = svd(data_combined_centered', 'econ');
    k_opt = 6; % Safe default or substitute with velicer_map output
    if k_opt > size(V,2), k_opt = size(V,2); end
    
    try
        X_dpca = cat(3, avgA, avgB); 
        X_dpca_centered = bsxfun(@minus, X_dpca, mean(X_dpca(:,:), 2));
        [W_dpca, ~, ~] = dpca(X_dpca_centered, k_opt); 
        W = W_dpca'; % Transpose to [k x Channels]
    catch 
        W = V(:, 1:k_opt)'; % Fallback to SVD
    end
    pc_labels = arrayfun(@(x) sprintf('dPC%d', x), 1:k_opt, 'UniformOutput', false);

    % 3. PROJECT TRIALS INTO dPC SPACE
    T_len = size(trialsA, 2);
    comp_trialsA = reshape(W * reshape(trialsA, num_ch, []), k_opt, T_len, nA);
    comp_trialsB = reshape(W * reshape(trialsB, num_ch, []), k_opt, T_len, nB);

    pool_trials_ch = cat(3, trialsA, trialsB);
    pool_trials_pc = cat(3, comp_trialsA, comp_trialsB);

    % 4. SLIDING WINDOW
    win_samples = round((window_size_ms / 1000) * fs);
    step_samples = round((step_size_ms / 1000) * fs);
    
    idx_start = find(time_ms_eeg >= t_win(1), 1, 'first');
    idx_end = find(time_ms_eeg <= t_win(2), 1, 'last');
    time_chunk = time_ms_eeg(idx_start:idx_end);
    
    start_idx = 1 : step_samples : (length(time_chunk) - win_samples + 1);
    num_windows = length(start_idx);

    sub_ch_A_all = nan(num_ch, num_ch, num_windows);
    sub_ch_B_all = nan(num_ch, num_ch, num_windows);
    sub_pc_A_all = nan(k_opt, k_opt, num_windows);
    sub_pc_B_all = nan(k_opt, k_opt, num_windows);
    window_centers = zeros(1, num_windows);

    for w = 1:num_windows
        w_s = idx_start + start_idx(w) - 1;
        w_e = w_s + win_samples - 1;
        window_centers(w) = mean(time_ms_eeg(w_s:w_e));

        % A. REAL CORRELATION (Trial-Averaged ERP)
        real_ch_A = corrcoef(mean(trialsA(:, w_s:w_e, :), 3, 'omitnan')');
        real_ch_B = corrcoef(mean(trialsB(:, w_s:w_e, :), 3, 'omitnan')');
        real_pc_A = corrcoef(mean(comp_trialsA(:, w_s:w_e, :), 3, 'omitnan')');
        real_pc_B = corrcoef(mean(comp_trialsB(:, w_s:w_e, :), 3, 'omitnan')');

        % B. BACKGROUND PERMUTATION (Fisher Z-Transformed)
        null_ch_A = nan(num_ch, num_ch, n_perms); null_ch_B = nan(num_ch, num_ch, n_perms);
        null_pc_A = nan(k_opt, k_opt, n_perms);   null_pc_B = nan(k_opt, k_opt, n_perms);

        pool_ch_slice = pool_trials_ch(:, w_s:w_e, :);
        pool_pc_slice = pool_trials_pc(:, w_s:w_e, :);

        for prm = 1:n_perms
            shuf_idx = randperm(total_trials);
            idxA = shuf_idx(1:nA); idxB = shuf_idx(nA+1:end);
            
            pseudo_ch_A = mean(pool_ch_slice(:,:,idxA), 3, 'omitnan');
            pseudo_ch_B = mean(pool_ch_slice(:,:,idxB), 3, 'omitnan');
            pseudo_pc_A = mean(pool_pc_slice(:,:,idxA), 3, 'omitnan');
            pseudo_pc_B = mean(pool_pc_slice(:,:,idxB), 3, 'omitnan');
            
            % Zero diagonal before atanh to avoid Infinity
            null_ch_A(:,:,prm) = atanh(corrcoef(pseudo_ch_A') .* ~eye(num_ch)); 
            null_ch_B(:,:,prm) = atanh(corrcoef(pseudo_ch_B') .* ~eye(num_ch));
            null_pc_A(:,:,prm) = atanh(corrcoef(pseudo_pc_A') .* ~eye(k_opt)); 
            null_pc_B(:,:,prm) = atanh(corrcoef(pseudo_pc_B') .* ~eye(k_opt));
        end

        % Inverse Fisher Transform
        bg_ch_A = tanh(mean(null_ch_A, 3, 'omitnan')); bg_ch_A(1:num_ch+1:end) = 1; 
        bg_ch_B = tanh(mean(null_ch_B, 3, 'omitnan')); bg_ch_B(1:num_ch+1:end) = 1;
        bg_pc_A = tanh(mean(null_pc_A, 3, 'omitnan')); bg_pc_A(1:k_opt+1:end) = 1; 
        bg_pc_B = tanh(mean(null_pc_B, 3, 'omitnan')); bg_pc_B(1:k_opt+1:end) = 1;

        % C. SUBTRACT
        sub_ch_A_all(:,:,w) = real_ch_A - bg_ch_A;
        sub_ch_B_all(:,:,w) = real_ch_B - bg_ch_B;
        sub_pc_A_all(:,:,w) = real_pc_A - bg_pc_A;
        sub_pc_B_all(:,:,w) = real_pc_B - bg_pc_B;
    end
end