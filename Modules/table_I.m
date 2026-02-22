pm = zeros(num_planes,1);
pa = zeros(num_planes,1);
b  = zeros(num_planes,1);

Ns = 2000;   % dense sampling for numerical integration

for i = 1:num_planes
    a1 = a_edges(i);
    a2 = a_edges(i+1);

    % Dense samples
    a_s = linspace(a1,a2,Ns);
    m_s = linspace(0,1,Ns);

    pm(i) = mean( sin(pi/2 * a_s) );

    [MM,AA] = meshgrid(m_s,a_s);
    dfa = (pi/2) * (1 + MM) .* cos(pi/2 * AA);
    pa(i) = mean(dfa(:));

    F_exact_sub = (1 + MM) .* sin(pi/2 * AA);
    F_plane_no_b = pm(i)*MM + pa(i)*AA;

    % MAE-minimizing bias (centers residual)
    b(i) = mean(F_exact_sub(:) - F_plane_no_b(:));
end

fprintf('\nTable I: Exact Plane Coefficients (Eq. 11 & 12)\n');
fprintf('-------------------------------------------------\n');
fprintf('Plane\tp_m\t\tp_a\t\tb\n');

for i = 1:num_planes
    fprintf('%d\t%.6f\t%.6f\t%.6f\n', ...
        i-1, pm(i), pa(i), b(i));
end
