function [encNet, decNet, priorNet, info] = trainEEG_iCAE(X_train_1D, C_train_1D, X_test_1D, C_test_1D, cfg)
    % X and C must be 4D tensors: [1, WindowLength, Channels, NumTrials]
    
    arguments
        X_train_1D double 
        C_train_1D double
        X_test_1D  double 
        C_test_1D  double
        cfg struct = struct()
    end
    
    % -------------------
    % Apply Default Configuration
    % -------------------
    if ~isfield(cfg, 'method'),         cfg.method = "icae"; end
    if ~isfield(cfg, 'bottleneckSize'), cfg.bottleneckSize = 8; end
    if ~isfield(cfg, 'epochs'),         cfg.epochs = 250; end
    if ~isfield(cfg, 'batchSize'),      cfg.batchSize = 32; end
    if ~isfield(cfg, 'learnRate'),      cfg.learnRate = 1e-4; end % Kept low for stability
    if ~isfield(cfg, 'patience'),       cfg.patience = 30; end
    if ~isfield(cfg, 'beta'),           cfg.beta = 0.1; end
    if ~isfield(cfg, 'warmupEpochs'),   cfg.warmupEpochs = 25; end
    
    WindowLength  = size(X_train_1D, 2);
    TotalChannels = size(X_train_1D, 3);
    NumClasses    = size(C_train_1D, 3);
    NumObs        = size(X_train_1D, 4);
    
    % ==========================================
    % 1. ENCODER NETWORK
    % ==========================================
    encInputDim = TotalChannels + NumClasses;
    
    lgraphEnc = layerGraph(imageInputLayer([1 WindowLength encInputDim], "Normalization", "none", "Name", "enc_in"));
    lgraphEnc = addLayers(lgraphEnc, [
        convolution2dLayer([1 15], 256, "Padding", "same", "Name", "enc_conv1")
        leakyReluLayer(0.01, "Name", "enc_leakyrelu1")
    ]);
    lgraphEnc = connectLayers(lgraphEnc, "enc_in", "enc_conv1");
    
    lgraphEnc = addLayers(lgraphEnc, convolution2dLayer([1 1], cfg.bottleneckSize, "Padding", "same", "Name", "enc_mu"));
    lgraphEnc = addLayers(lgraphEnc, convolution2dLayer([1 1], cfg.bottleneckSize, "Padding", "same", "Name", "enc_logvar"));
    lgraphEnc = connectLayers(lgraphEnc, "enc_leakyrelu1", "enc_mu");
    lgraphEnc = connectLayers(lgraphEnc, "enc_leakyrelu1", "enc_logvar");
    
    encNet = dlnetwork(lgraphEnc, dlarray(zeros(1, WindowLength, encInputDim, 1), 'SSCB'));
    
    % ==========================================
    % 2. DECODER NETWORK (Removed tanh layer to allow normal EEG scales)
    % ==========================================
    decInputDim = cfg.bottleneckSize + NumClasses;
    
    lgraphDec = layerGraph(imageInputLayer([1 WindowLength decInputDim], "Normalization", "none", "Name", "dec_in"));
    lgraphDec = addLayers(lgraphDec, [
        convolution2dLayer([1 15], 256, "Padding", "same", "Name", "dec_conv1")
        leakyReluLayer(0.01, "Name", "dec_leakyrelu1")
        convolution2dLayer([1 1], TotalChannels, "Padding", "same", "Name", "reconstruction")
    ]);
    lgraphDec = connectLayers(lgraphDec, "dec_in", "dec_conv1");
    
    decNet = dlnetwork(lgraphDec, dlarray(zeros(1, WindowLength, decInputDim, 1), 'SSCB'));
    
    % ==========================================
    % 3. PRIOR NETWORK
    % ==========================================
    priorNet = [];
    if cfg.method == "icae"
        lgraphPrior = layerGraph(imageInputLayer([1 WindowLength NumClasses], "Normalization", "none", "Name", "prior_in"));
        lgraphPrior = addLayers(lgraphPrior, [
            convolution2dLayer([1 1], 16, "Padding", "same", "Name", "prior_conv1")
            reluLayer("Name", "prior_relu1")
        ]);
        lgraphPrior = connectLayers(lgraphPrior, "prior_in", "prior_conv1");
        
        lgraphPrior = addLayers(lgraphPrior, convolution2dLayer([1 1], cfg.bottleneckSize, "Padding", "same", "Name", "prior_mu"));
        lgraphPrior = addLayers(lgraphPrior, convolution2dLayer([1 1], cfg.bottleneckSize, "Padding", "same", "Name", "prior_logvar"));
        lgraphPrior = connectLayers(lgraphPrior, "prior_relu1", "prior_mu");
        lgraphPrior = connectLayers(lgraphPrior, "prior_relu1", "prior_logvar");
        
        priorNet = dlnetwork(lgraphPrior, dlarray(zeros(1, WindowLength, NumClasses, 1), 'SSCB'));
    end

    % ==========================================
    % 4. TRAINING LOOP SETUP
    % ==========================================
    X_dl = dlarray(X_train_1D, 'SSCB'); 
    C_dl = dlarray(C_train_1D, 'SSCB');
    X_test_dl = dlarray(X_test_1D, 'SSCB');
    C_test_dl = dlarray(C_test_1D, 'SSCB');
    
    trailingAvgEnc = []; trailingAvgSqEnc = [];
    trailingAvgDec = []; trailingAvgSqDec = [];
    trailingAvgPrior = []; trailingAvgSqPrior = [];
    
    numBatches = floor(NumObs / cfg.batchSize);
    
    info.lossHistory = zeros(cfg.epochs, 1);
    info.valLossHistory = zeros(cfg.epochs, 1);
    
    bestValLoss = inf;
    patienceCounter = 0;
    
    info.bestEncNet = encNet;
    info.bestDecNet = decNet;
    info.bestPriorNet = priorNet;
    
    % ==========================================
    % 5. EXECUTE TRAINING
    % ==========================================
    for epoch = 1:cfg.epochs
        epochLoss = 0;
        
        % KL Annealing (Dynamic Beta)
        if cfg.warmupEpochs > 1
            fraction = min(1, (epoch - 1) / (cfg.warmupEpochs - 1));
            current_beta = cfg.beta * fraction;
        else
            current_beta = cfg.beta;
        end
        
        % Shuffle Data
        idx = randperm(NumObs);
        X_shuffled = X_dl(:, :, :, idx);
        C_shuffled = C_dl(:, :, :, idx);
        
        for batch = 1:numBatches
            idxBatch = (batch-1)*cfg.batchSize + 1 : batch*cfg.batchSize;
            XBatch = X_shuffled(:, :, :, idxBatch);
            CBatch = C_shuffled(:, :, :, idxBatch);
            
            [gradEnc, gradDec, gradPrior, loss] = dlfeval(@modelLoss, encNet, decNet, priorNet, XBatch, CBatch, cfg.method, current_beta);
            
            % --- DEFENSE #1: GRADIENT CLIPPING ---
            % This strictly limits gradients to between -1 and 1, preventing the weights from exploding
            gradThreshold = 1.0;
            gradEnc = dlupdate(@(g) min(max(g, -gradThreshold), gradThreshold), gradEnc);
            gradDec = dlupdate(@(g) min(max(g, -gradThreshold), gradThreshold), gradDec);
            if cfg.method == "icae"
                gradPrior = dlupdate(@(g) min(max(g, -gradThreshold), gradThreshold), gradPrior);
            end
            % -------------------------------------
            
            [encNet, trailingAvgEnc, trailingAvgSqEnc] = adamupdate(encNet, gradEnc, trailingAvgEnc, trailingAvgSqEnc, epoch, cfg.learnRate);
            [decNet, trailingAvgDec, trailingAvgSqDec] = adamupdate(decNet, gradDec, trailingAvgDec, trailingAvgSqDec, epoch, cfg.learnRate);
            
            if cfg.method == "icae"
                [priorNet, trailingAvgPrior, trailingAvgSqPrior] = adamupdate(priorNet, gradPrior, trailingAvgPrior, trailingAvgSqPrior, epoch, cfg.learnRate);
            end
            
            epochLoss = epochLoss + extractdata(loss);
        end
        
        % Record and Validate
        info.lossHistory(epoch) = epochLoss / numBatches;
        valLoss = computeValidationLoss(encNet, decNet, priorNet, X_test_dl, C_test_dl, cfg.method, current_beta);
        info.valLossHistory(epoch) = extractdata(valLoss);
        
        fprintf("Epoch %d/%d (Beta: %.4f) - Train Loss: %.4f | Val Loss: %.4f\n", ...
                epoch, cfg.epochs, current_beta, info.lossHistory(epoch), info.valLossHistory(epoch));
                
        % Early Stopping Logic
        currentValLoss = info.valLossHistory(epoch);
        
        if isnan(currentValLoss) || isinf(currentValLoss)
            warning('Validation loss hit NaN/Inf. Gradient clipping actively suppressing explosion.');
            patienceCounter = patienceCounter + 1;
        elseif currentValLoss < bestValLoss
            bestValLoss = currentValLoss;
            patienceCounter = 0;
            info.bestEncNet = encNet;
            info.bestDecNet = decNet;
            info.bestPriorNet = priorNet;
        else
            patienceCounter = patienceCounter + 1;
        end
        
        if patienceCounter >= cfg.patience
            fprintf("Early stopping triggered! Validation loss hasn't improved for %d epochs.\n", cfg.patience);
            info.lossHistory = info.lossHistory(1:epoch);
            info.valLossHistory = info.valLossHistory(1:epoch);
            
            encNet = info.bestEncNet; 
            decNet = info.bestDecNet; 
            priorNet = info.bestPriorNet;
            break;
        end
    end
end

% -------------------
% Helper: Model Loss Function
% -------------------
function [gradEnc, gradDec, gradPrior, loss] = modelLoss(encNet, decNet, priorNet, X, C, method, beta)
    XC = cat(3, X, C);
    [mu_q, logvar_q] = forward(encNet, XC);
    
    if method == "icae"
        [mu_p, logvar_p] = forward(priorNet, C);
    else
        mu_p = zeros(size(mu_q), 'like', mu_q);
        logvar_p = zeros(size(logvar_q), 'like', logvar_q);
    end
    
    % --- DEFENSE #2: LOG-VAR CLAMPING ---
    % Force the log variance to stay inside a numerically stable range (-10 to 10)
    % This stops exp(logvar) from turning into Infinity or Zero!
    logvar_q = max(min(logvar_q, 10), -10);
    logvar_p = max(min(logvar_p, 10), -10);
    % ------------------------------------
    
    epsilon = randn(size(mu_q), 'like', mu_q);
    Z = mu_q + exp(0.5 * logvar_q) .* epsilon;
    
    ZC = cat(3, Z, C);
    X_hat = forward(decNet, ZC);
    
    mseLoss = mean((X_hat - X).^2, 'all'); 
    var_q = exp(logvar_q);
    var_p = exp(logvar_p);
    
    klDivergence = 0.5 * sum(logvar_p - logvar_q + (var_q + (mu_q - mu_p).^2)./var_p - 1, 3);
    klLoss = mean(klDivergence, 'all');
    
    loss = mseLoss + beta * klLoss;
    
    if method == "icae"
        [gradEnc, gradDec, gradPrior] = dlgradient(loss, encNet.Learnables, decNet.Learnables, priorNet.Learnables);
    else
        [gradEnc, gradDec] = dlgradient(loss, encNet.Learnables, decNet.Learnables);
        gradPrior = [];
    end
end

% -------------------
% Helper: Validation Loss
% -------------------
function valLoss = computeValidationLoss(encNet, decNet, priorNet, X, C, method, beta)
    XC = cat(3, X, C);
    [mu_q, logvar_q] = forward(encNet, XC);
    
    if method == "icae"
        [mu_p, logvar_p] = forward(priorNet, C);
    else
        mu_p = zeros(size(mu_q), 'like', mu_q);
        logvar_p = zeros(size(logvar_q), 'like', logvar_q);
    end
    
    % Apply the same Log-Var clamping to validation to prevent NaN checks
    logvar_q = max(min(logvar_q, 10), -10);
    logvar_p = max(min(logvar_p, 10), -10);
    
    Z = mu_q; 
    ZC = cat(3, Z, C);
    X_hat = forward(decNet, ZC);
    
    mseLoss = mean((X_hat - X).^2, 'all');
    var_q = exp(logvar_q);
    var_p = exp(logvar_p);
    
    klDivergence = 0.5 * sum(logvar_p - logvar_q + (var_q + (mu_q - mu_p).^2)./var_p - 1, 3);
    klLoss = mean(klDivergence, 'all');
    
    valLoss = mseLoss + beta * klLoss;
end