clear; clc; close all;

F = @(m,a) (1+m).*2.^(-a);

m_vals = linspace(0,1,256);
a_vals = linspace(0,1,256);
[M,A] = meshgrid(m_vals,a_vals);

F_exact = F(M,A);
Y = F_exact(:);

X_base = [ones(numel(M),1) M(:) A(:)];
coef_base = X_base \ Y;

F_base = reshape(X_base * coef_base, size(M));

bit_list = [0 2 4 6 8];   % 0 = baseline

results_MRED = zeros(size(bit_list));
results_NMED = zeros(size(bit_list));

for k = 1:length(bit_list)

    bits = bit_list(k);

    if bits == 0
        % Baseline
        F_model = F_base;
    else
        % Quantize angle ONLY for cross term
        A_q = round(A*(2^bits - 1))/(2^bits - 1);

        X_prop = [ones(numel(M),1) ...
                  M(:) ...
                  A(:) ...
                  (M(:).*A_q(:))];

        coef_prop = X_prop \ Y;

        F_model = reshape(X_prop * coef_prop, size(M));
    end

    err = F_model - F_exact;

    % MRED
    MRED = mean(abs(err(:) ./ F_exact(:)));

    % NMED
    NMED = sqrt(mean(err(:).^2)) / max(F_exact(:));

    results_MRED(k) = MRED * 100;
    results_NMED(k) = NMED * 100;

end

fprintf('\n===== ERROR TABLE (Error / ×10^-2) =====\n\n');

for k = 1:length(bit_list)

    if bit_list(k) == 0
        label = 'Baseline';
    else
        label = sprintf('Proposed_%d', bit_list(k));
    end

    fprintf('%s:\n', label);
    fprintf('MRED = %.4f\n', results_MRED(k));
    fprintf('NMED = %.4f\n\n', results_NMED(k));
end