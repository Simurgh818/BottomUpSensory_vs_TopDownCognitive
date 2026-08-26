function [num_factors, map_values] = velicer_map(data, quiet)
% VELICER_MAP Determines the number of factors to retain using Velicer's MAP test.
%
% Inputs:
%   data  - An (n x p) matrix of raw data (n observations, p variables).
%   quiet - Boolean flag. Set to true to suppress plotting (default: false).

    if nargin < 2
        quiet = false; % Default is to show the plot
    end

    % 1. Compute the correlation matrix
    R = corrcoef(data);
    p = size(R, 1);
    
    % 2. Perform Eigenvalue Decomposition
    [eigvec, eigval_mat] = eig(R);
    eigval = diag(eigval_mat);
    
    % Sort eigenvalues and eigenvectors in descending order
    [eigval, idx] = sort(eigval, 'descend');
    eigvec = eigvec(:, idx);
    
    % Initialize array to store the MAP values
    map_values = zeros(1, p);
    
    % 3. Step 0: Calculate average squared correlation (no factors partialled out)
    R_off_diag = R - eye(p);
    map_values(1) = sum(R_off_diag(:).^2) / (p * (p - 1));
    
    % 4. Loop to partial out 1 to p-1 components
    for m = 1:(p-1)
        A = eigvec(:, 1:m) * diag(sqrt(eigval(1:m)));
        part_cov = R - (A * A');
        d = diag(part_cov);
        inv_sqrt_d = diag(1 ./ sqrt(d));
        R_partial = inv_sqrt_d * part_cov * inv_sqrt_d;
        R_partial_off = R_partial - eye(p);
        map_values(m+1) = sum(R_partial_off(:).^2) / (p * (p - 1));
    end
    
    % 5. Find the minimum MAP value
    [min_val, min_idx] = min(map_values);
    num_factors = min_idx - 1;
    
    % 6. Optional: Plot the results to visualize the minimum
    if ~quiet
        figure('Position',[50 50 800 500]);
        plot(0:(p-1), map_values, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
        hold on;
        plot(num_factors, map_values(min_idx), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r'); 
        coordinate_text = sprintf('  Min: (%d, %.4f)', num_factors, min_val);
        text(num_factors, 0.2, coordinate_text, 'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'center', 'FontSize', 26, 'FontWeight', 'bold', 'Color', 'r');
        xlabel('Number of Factors Retained');
        ylabel('Avg. Sqrd. Partial Corr.');
        title('Velicer''s MAP Test');
        legend('MAP Values', 'Minimum (Optimal Factors)');
        set(gca,'FontSize',24);
        grid on; hold off;
    end
end
