clear; clc;

% Grid
m = linspace(0,1,256);
a = linspace(0,1,256);
[M,A] = meshgrid(m,a);

% True CM surface
F_exact = (1 + M).*2.^(-A);

X_base = [M(:) A(:) ones(numel(M),1)];
y = F_exact(:);

coeff_base = X_base \ y;

F_base = reshape(X_base*coeff_base, size(F_exact));

X_prop = [M(:) A(:) (M(:).*A(:)) ones(numel(M),1)];
coeff_prop = X_prop \ y;

F_prop = reshape(X_prop*coeff_prop, size(F_exact));
eps_val = 1e-12;

RED_base = abs((F_base - F_exact)./(F_exact + eps_val));
RED_prop = abs((F_prop - F_exact)./(F_exact + eps_val));

MRED_base = mean(RED_base(:));
MRED_prop = mean(RED_prop(:));

NMED_base = mean(abs(F_base(:) - F_exact(:))) / mean(F_exact(:));
NMED_prop = mean(abs(F_prop(:) - F_exact(:))) / mean(F_exact(:));

fprintf('\n===== PAPER-STYLE ERROR METRICS (×10^-2) =====\n\n');

fprintf('Baseline CM:\n');
fprintf('MRED = %.4f\n', MRED_base*100);
fprintf('NMED = %.4f\n\n', NMED_base*100);

fprintf('Proposed CM (Interaction + cross-term):\n');
fprintf('MRED = %.4f\n', MRED_prop*100);
fprintf('NMED = %.4f\n\n', NMED_prop*100);

fprintf('Improvement:\n');
fprintf('MRED reduction = %.2f %%\n', ...
    100*(MRED_base - MRED_prop)/MRED_base);
fprintf('NMED reduction = %.2f %%\n', ...
    100*(NMED_base - NMED_prop)/NMED_base);
