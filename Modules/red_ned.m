clear; clc; close all;

Nm = 300;
Na = 300;

m = linspace(0,1,Nm);
a = linspace(0,1,Na);

[M,A] = meshgrid(m,a);


F_exact = (1 + M) .* sin(pi/2 * A);

pm  = 0.95;
pa  = 0.78;
pma = 0.32;
b   = 0.01;

F_approx = pm*M + pa*A + pma*(M.*A) + b;

Error = F_approx - F_exact;

RED = Error ./ (F_exact + eps);
NED = Error ./ max(abs(F_exact(:)));

threshold = 0.05 * max(F_exact(:));
RED_masked = RED;
RED_masked(abs(F_exact) < threshold) = NaN;

figure;
imagesc(a, m, RED_masked.');
set(gca,'YDir','normal');
colormap jet;
colorbar;

xlabel('Angle (a)');
ylabel('Mantissa (m)');
title('(a) Relative Error Distribution (RED)');

hold on;
contour(a, m, RED_masked.', 10, 'k');
hold off;

figure;
imagesc(a, m, NED.');
set(gca,'YDir','normal');
colormap jet;
colorbar;

xlabel('Angle (a)');
ylabel('Mantissa (m)');
title('(b) Normalized Error Distribution (NED)');

hold on;
contour(a, m, NED.', 10, 'k');
hold off;
