clear; clc; close all;

Nm = 256;
Na = 256;

m = linspace(0,1,Nm);
a = linspace(0,1,Na);

[M,A] = meshgrid(m,a);

F_exact = (1 + M) .* sin(pi/2 * A);

num_planes = 8;
a_edges = linspace(0,1,num_planes+1);

fprintf('\nUniform plane boundaries (Stage 4):\n');
for i = 1:num_planes
    fprintf('Plane %d: a ∈ [%.4f , %.4f)\n', ...
        i-1, a_edges(i), a_edges(i+1));
end

pm  = zeros(num_planes,1);
pa  = zeros(num_planes,1);
pma = zeros(num_planes,1);
b   = zeros(num_planes,1);

Ns = 120;        % samples per dimension
max_iter = 15;
eps = 1e-6;

for i = 1:num_planes
    a1 = a_edges(i);
    a2 = a_edges(i+1);

    % Sample points inside plane
    m_s = linspace(0,1,Ns);
    a_s = linspace(a1,a2,Ns);
    [MM,AA] = meshgrid(m_s,a_s);

    % Exact surface
    F = (1 + MM) .* sin(pi/2 * AA);

    % Interaction-aware regression matrix
    % f ≈ pm*m + pa*a + pma*(m*a) + b
    X = [ MM(:), ...
          AA(:), ...
          (MM(:).*AA(:)), ...
          ones(numel(MM),1) ];

    y = F(:);

    % Initial least-squares solution
    theta = X \ y;

    % IRLS loop (MAE-oriented)
    for k = 1:max_iter
        r = y - X*theta;
        w = 1 ./ (abs(r) + eps);

        Xw = X .* w;
        theta = (X' * Xw) \ (Xw' * y);
    end

    pm(i)  = theta(1);
    pa(i)  = theta(2);
    pma(i) = theta(3);
    b(i)   = theta(4);
end

F_approx = zeros(size(F_exact));

for i = 1:num_planes
    mask = (A >= a_edges(i)) & (A < a_edges(i+1));
    F_approx(mask) = ...
        pm(i)*M(mask) + ...
        pa(i)*A(mask) + ...
        pma(i)*(M(mask).*A(mask)) + ...
        b(i);
end

abs_err = abs(F_exact - F_approx);

MAE  = mean(abs_err(:));
MRED = mean(abs_err(:)) / mean(abs(F_exact(:)));

fprintf('\nStage 4 (Interaction-Aware) Results:\n');
fprintf('MAE  = %.4e\n', MAE);
fprintf('MRED = %.4e\n', MRED);

figure;
surf(M, A, F_exact, 'EdgeColor','none');
title('Exact Surface f(m,a)');
xlabel('m'); ylabel('a'); zlabel('f');

figure;
surf(M, A, F_approx, 'EdgeColor','none');
title('Interaction-Aware Approximate Surface');
xlabel('m'); ylabel('a'); zlabel('f');

hold on;

for i = 2:num_planes
    a_b = a_edges(i);
    m_line = linspace(0,1,200);
    a_line = a_b * ones(size(m_line));

    % Evaluate approximation exactly on the boundary
    % Find which plane index applies just below the boundary
    idx = i-1;

    f_line = ...
        pm(idx)*m_line + ...
        pa(idx)*a_line + ...
        pma(idx)*(m_line.*a_line) + ...
        b(idx);

    plot3(m_line, a_line, f_line, ...
          'k', 'LineWidth', 2);
end

hold off;


figure;
surf(M, A, abs_err, 'EdgeColor','none');
title('Absolute Error Surface');
xlabel('m'); ylabel('a'); zlabel('|Error|');

figure;
histogram(abs_err(:),100);
title('Error Distribution');
xlabel('Absolute Error'); ylabel('Count');

fprintf('\nInteraction-Aware Coefficients\n');
fprintf('Plane\tpm\t\tpa\t\tpma\t\tb\n');

for i = 1:num_planes
    fprintf('%d\t%.6f\t%.6f\t%.6f\t%.6f\n', ...
        i-1, pm(i), pa(i), pma(i), b(i));
end

