clear; clc;

B = 4;   % <-- CHANGE THIS to any of 2,4,6,8 bits

F = @(m,a) (1+m).*2.^(-a);

m_vals = linspace(0,1,256);
a_vals = linspace(0,1,256);
[M,A] = meshgrid(m_vals,a_vals);
F_exact = F(M,A);
edges = linspace(0,1,9);

RMSE_base = zeros(8,1);
RMSE_int  = zeros(8,1);
Improvement = zeros(8,1);

if B > 0
    A_q_full = round(A*(2^B - 1))/(2^B - 1);
else
    A_q_full = zeros(size(A));   % no interaction
end


for p = 1:8

    mask = (A >= edges(p)) & (A < edges(p+1));

    m_plane = M(mask);
    a_plane = A(mask);
    a_q_plane = A_q_full(mask);
    y_plane = F_exact(mask);

    %% Baseline

    X_base = [m_plane(:) a_plane(:) ones(length(m_plane),1)];
    coef_base = X_base \ y_plane(:);
    y_hat_base = X_base * coef_base;

    RMSE_base(p) = sqrt(mean((y_hat_base - y_plane(:)).^2));

    %% Interaction-aware

    if B > 0
        X_int = [m_plane(:) a_plane(:) ...
                 (m_plane(:).*a_q_plane(:)) ...
                 ones(length(m_plane),1)];
    else
        X_int = X_base;   % same as baseline
    end

    coef_int = X_int \ y_plane(:);
    y_hat_int = X_int * coef_int;

    RMSE_int(p) = sqrt(mean((y_hat_int - y_plane(:)).^2));

    Improvement(p) = ...
        (RMSE_base(p) - RMSE_int(p)) ...
        / RMSE_base(p) * 100;

end

fprintf('\n===== Plane-wise RMSE (B = %d bits) =====\n\n',B);

for p = 1:8
    fprintf('Plane %d:\n',p-1);
    fprintf('RMSE (No Interaction)   = %.6e\n',RMSE_base(p));
    fprintf('RMSE (With Interaction) = %.6e\n',RMSE_int(p));
    fprintf('Improvement = %.2f %%\n\n',Improvement(p));
end

fprintf('Average Improvement = %.2f %%\n',mean(Improvement));