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

bit_widths = [6, 8];

MAE_fp   = zeros(length(bit_widths),1);
MRED_fp = zeros(length(bit_widths),1);

for bw_idx = 1:length(bit_widths)

    B = bit_widths(bw_idx);
    Q = 2^B;

    % Quantize angle 'a' to B bits
    A_q = round(A * Q) / Q;

    % Apply interaction-aware plane approximation
    F_approx_q = zeros(size(F_exact));

    for i = 1:num_planes
        mask = (A_q >= a_edges(i)) & (A_q < a_edges(i+1));

        F_approx_q(mask) = ...
            pm(i)*M(mask) + ...
            pa(i)*A_q(mask) + ...
            pma(i)*(M(mask).*A_q(mask)) + ...
            b(i);
    end

    % Error metrics
    abs_err = abs(F_exact - F_approx_q);

    MAE_fp(bw_idx)  = mean(abs_err(:));
    MRED_fp(bw_idx) = mean(abs_err(:)) / mean(abs(F_exact(:)));

    fprintf('\nStage 4B Results (Angle = %d bits):\n', B);
    fprintf('MAE  = %.4e\n', MAE_fp(bw_idx));
    fprintf('MRED = %.4e\n', MRED_fp(bw_idx));
end


B = 8;
Q = 2^B;

% Quantize angle to 8 bits
A_q = round(A * Q) / Q;

% Compute bit-slice m·a
MA_bitslice = bit_slice_mult(M, A_q, B);

% Apply interaction-aware plane approximation (bit-slice)
F_approx_bs = zeros(size(F_exact));

for i = 1:num_planes
    mask = (A_q >= a_edges(i)) & (A_q < a_edges(i+1));

    F_approx_bs(mask) = ...
        pm(i)*M(mask) + ...
        pa(i)*A_q(mask) + ...
        pma(i)*MA_bitslice(mask) + ...
        b(i);
end

% Error metrics
abs_err_bs = abs(F_exact - F_approx_bs);

MAE_bs  = mean(abs_err_bs(:));
MRED_bs = mean(abs_err_bs(:)) / mean(abs(F_exact(:)));

fprintf('\nStage 4C Results (8-bit, bit-slice m·a):\n');
fprintf('MAE  = %.4e\n', MAE_bs);
fprintf('MRED = %.4e\n', MRED_bs);


function ma = bit_slice_mult(m, a_q, B)
% Bit-slice multiplication: m * a_q
% m   : mantissa (same size as a_q)
% a_q : fixed-point angle in [0,1)
% B   : number of fractional bits (here B = 8)

    ma = zeros(size(m));

    for k = 1:B
        % Extract k-th fractional bit of a_q
        a_bit = bitget(floor(a_q * 2^B), B - k + 1);

        % Accumulate shifted m where bit is 1
        ma = ma + a_bit .* (m / 2^k);
    end
end
