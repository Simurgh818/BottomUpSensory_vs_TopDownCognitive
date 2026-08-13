function [sub_ch_A_all, sub_ch_B_all, sub_pc_A_all, sub_pc_B_all, window_centers, k_opt, pc_labels] = ...
    compute_dynamic_connectivity(trialsA_raw, trialsB_raw, t_win, time_ms_eeg, fs, window_size_ms, step_size_ms, bg_windows)
    
    nA = size(trialsA_raw, 3);
    nB = size(trialsB_raw, 3);
    num_ch = size(trialsA_raw, 1);
    
    % --- 1. CALCULATE RESTING STATE BACKGROUND (Trial-Averaged) ---
    % First, get the trial-averaged ERP for each condition
    erp_A = mean(trialsA_raw, 3, 'omitnan');
    erp_B = mean(trialsB_raw, 3, 'omitnan');
    
    % Concatenate the start and end resting states across the ERP
    concat_bg_A_ch = [];
    concat_bg_B_ch = [];
    for w = 1:size(bg_windows, 1)
        idx_bg = find(time_ms_eeg >= bg_windows(w, 1) & time_ms_eeg <= bg_windows(w, 2));
        concat_bg_A_ch = [concat_bg_A_ch, erp_A(:, idx_bg)];
        concat_bg_B_ch = [concat_bg_B_ch, erp_B(:, idx_bg)];
    end
    
    % Baseline Channel correlation matrices
    bg_ch_A = corrcoef(concat_bg_A_ch');
    bg_ch_B = corrcoef(concat_bg_B_ch');
    
    % Baseline PC correlation matrices
    bg_pc_A = corrcoef((W * concat_bg_A_ch)');
    bg_pc_B = corrcoef((W * concat_bg_B_ch)');
    
    % --- 2. DIMENSIONALITY REDUCTION (SHARED MANIFOLD via dPCA) ---
    avgA = mean(trialsA_raw, 3, 'omitnan');
    avgB = mean(trialsB_raw, 3, 'omitnan');
    data_combined = [avgA, avgB];
    data_combined_centered = data_combined - mean(data_combined, 2);
    
    [~, ~, V] = svd(data_combined_centered', 'econ');
    k_opt = 6; 
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
    
    % Baseline PC correlation matrices
    bg_pc_A = corrcoef((W * concat_bg_A)');
    bg_pc_B = corrcoef((W * concat_bg_B)');
    
    % --- 3. SLIDING WINDOW EVALUATION ---
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
        
        % A. Calculate REAL correlation for this window on the Trial Average (ERP)
        act_A = mean(trialsA_raw(:, w_s:w_e, :), 3, 'omitnan');
        act_B = mean(trialsB_raw(:, w_s:w_e, :), 3, 'omitnan');
        
        real_ch_A = corrcoef(act_A');
        real_ch_B = corrcoef(act_B');
        
        real_pc_A = corrcoef((W * act_A)');
        real_pc_B = corrcoef((W * act_B)');
        
        % Subtract Stabilized Resting Background
        sub_ch_A_all(:,:,w) = real_ch_A - bg_ch_A;
        sub_ch_B_all(:,:,w) = real_ch_B - bg_ch_B;
        sub_pc_A_all(:,:,w) = real_pc_A - bg_pc_A;
        sub_pc_B_all(:,:,w) = real_pc_B - bg_pc_B;
    end
end