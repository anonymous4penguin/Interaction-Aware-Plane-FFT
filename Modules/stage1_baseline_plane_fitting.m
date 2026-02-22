% Uniform 8-plane partitioning
% Analytical coefficients
% Reference paper reproduction

clear; clc; close all;

Nm = 256;   % mantissa resolution
Na = 256;   % angle resolution

m = linspace(0,1,Nm);
a = linspace(0,1,Na);

[M,A] = meshgrid(m,a);

F_exact = (1 + M) .* sin(pi/2 * A);

num_planes = 8;
a_edges = linspace(0,1,num_planes+1);

pm = zeros(num_planes,1);
pa = zeros(num_planes,1);
b  = zeros(num_planes,1);

for i = 1:num_planes
    a1 = a_edges(i);
    a2 = a_edges(i+1);

    a_s = linspace(a1,a2,200);
    m_s = linspace(0,1,200);

    % p_m average slope in m-direction
    pm(i) = mean( sin(pi/2 * a_s) );

    % p_a average slope in a-direction
    [MM,AA] = meshgrid(m_s,a_s);
    dfa = (pi/2) * (1 + MM) .* cos(pi/2 * AA);
    pa(i) = mean(dfa(:));

    % b MAE-minimizing offset
    F_sub = (1 + MM) .* sin(pi/2 * AA);
    F_lin = pm(i)*MM + pa(i)*AA;
    b(i) = mean(F_sub(:) - F_lin(:));
end

F_approx = zeros(size(F_exact));

for i = 1:num_planes
    mask = (A >= a_edges(i)) & (A < a_edges(i+1));
    F_approx(mask) = pm(i)*M(mask) + pa(i)*A(mask) + b(i);
end

abs_err = abs(F_exact - F_approx);

MAE  = mean(abs_err(:));
MRED = mean(abs_err(:)) / mean(abs(F_exact(:)));


fprintf("Baseline MAE  = %.4e\n", MAE);
fprintf("Baseline MRED = %.4e\n", MRED);

figure;
surf(M, A, F_exact, 'EdgeColor', 'none');
title('Exact Surface f(m,a)');
xlabel('m'); ylabel('a'); zlabel('f');

figure;
surf(M, A, F_approx, 'EdgeColor', 'none');
title('Approximate Surface (8 Uniform Planes)');
xlabel('m'); ylabel('a'); zlabel('f');

figure;
surf(M, A, abs_err, 'EdgeColor', 'none');
title('Absolute Error Surface');
xlabel('m'); ylabel('a'); zlabel('|Error|');

figure;
histogram(abs_err(:), 100);
title('Error Distribution');
xlabel('Absolute Error'); ylabel('Count');

fprintf('\nTable-I Exact Plane Coefficients:\n');
fprintf('Plane\tpm\t\tpa\t\tb\n');

for i = 1:num_planes
    fprintf('%d\t%.6f\t%.6f\t%.6f\n', ...
        i-1, pm(i), pa(i), b(i));
end


%Non-uniform partitioning alone has NOT guarantee lower MAE when coefficients are 
% computed using slope averaging.

%Non-uniform partitioning requires coefficient optimization to be effective.
%This proves coefficient computation and partitioning are coupled