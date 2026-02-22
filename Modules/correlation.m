clear; clc; close all;

% Grid
Nm = 300;
Na = 300;

m = linspace(0,1,Nm);
a = linspace(0,1,Na);

[M,A] = meshgrid(m,a);

% Exact interaction component
Interaction = M .* sin(pi/2 * A);

% Normalize for visualization
Interaction_norm = Interaction / max(Interaction(:));

figure;
imagesc(m, a, Interaction_norm);
set(gca,'YDir','normal');
colormap jet;
colorbar;

xlabel('Mantissa (m)');
ylabel('Angle (a)');
title('Heatmap of Interaction Term  m · sin(\pi/2 · a)');

hold on;
contour(m, a, Interaction_norm, 10, 'k');
hold off;
