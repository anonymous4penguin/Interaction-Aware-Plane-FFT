clear; clc; close all;

Nm = 256;
Na = 256;

m = linspace(0,1,Nm);
a = linspace(0,1,Na);

[M,A] = meshgrid(m,a);

F_exact = (1 + M) .* sin(pi/2 * A);
num_planes = 8;

Ns = 120;
a_dense = linspace(0,1,Ns);

curv = sin(pi/2 * a_dense);                % curvature proxy
cum_curv = cumtrapz(a_dense, curv);
cum_curv = cum_curv / cum_curv(end);       % normalize

a_edges = zeros(num_planes+1,1);
a_edges(1) = 0;
a_edges(end) = 1;

for i = 2:num_planes
    target = (i-1)/num_planes;
    [~,idx] = min(abs(cum_curv - target));
    a_edges(i) = a_dense(idx);
end

fprintf('\nStage 2: Non-uniform plane boundaries\n');
for i = 1:num_planes
    fprintf('Plane %d: a ∈ [%.4f , %.4f)\n', ...
        i-1, a_edges(i), a_edges(i+1));
end

pm = zeros(num_planes,1);
pa = zeros(num_planes,1);
b  = zeros(num_planes,1);

Ns = 200;        % samples per dimension
max_iter = 15;
eps = 1e-6;

for i = 1:num_planes
    a1 = a_edges(i);
    a2 = a_edges(i+1);

    % Sample points in plane
    m_s = linspace(0,1,Ns);
    a_s = linspace(a1,a2,Ns);
    [MM,AA] = meshgrid(m_s,a_s);

    F = (1 + MM) .* sin(pi/2 * AA);

    X = [MM(:), AA(:), ones(numel(MM),1)];
    y = F(:);

    % Initial LS solution
    theta = X \ y;

    for k = 1:max_iter
        r = y - X*theta;
        w = 1 ./ (abs(r) + eps);   % weights

        % Weighted least squares without forming W
        Xw = X .* w;               % each row weighted
        theta = (X' * Xw) \ (Xw' * y);
    end

    pm(i) = theta(1);
    pa(i) = theta(2);
    b(i)  = theta(3);
end

F_approx = zeros(size(F_exact));

for i = 1:num_planes
    mask = (A >= a_edges(i)) & (A < a_edges(i+1));
    F_approx(mask) = pm(i)*M(mask) + pa(i)*A(mask) + b(i);
end

abs_err = abs(F_exact - F_approx);

MAE  = mean(abs_err(:));
MRED = mean(abs_err(:)) / mean(abs(F_exact(:)));

fprintf('\nStage 3B_B Results:\n');
fprintf('MAE  = %.4e\n', MAE);
fprintf('MRED = %.4e\n', MRED);

figure;
surf(M, A, F_exact, 'EdgeColor','none');
title('Exact Surface f(m,a)');
xlabel('m'); ylabel('a'); zlabel('f');

figure;
surf(M, A, F_approx, 'EdgeColor','none');
title('Approximate Surface (8 Non-Uniform Planes)');
xlabel('m'); ylabel('a'); zlabel('f');

figure;
surf(M, A, abs_err, 'EdgeColor','none');
title('Absolute Error Surface');
xlabel('m'); ylabel('a'); zlabel('|Error|');

figure;
histogram(abs_err(:),100);
title('Error Distribution');
xlabel('Absolute Error'); ylabel('Count');

fprintf('\nTable-I Coefficients (Stage 2 – Non-uniform)\n');
fprintf('Plane\tpm\t\tpa\t\tb\n');

for i = 1:num_planes
    fprintf('%d\t%.6f\t%.6f\t%.6f\n', ...
        i-1, pm(i), pa(i), b(i));
end
