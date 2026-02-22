clear; clc; rng(1);

N = 1024;

%% Load EEG
eeg_all = [];
for k = 0:35
    data = readmatrix(sprintf('s%02d.csv',k));
    eeg_all = [eeg_all; data(:,1)];
end
eeg_all = eeg_all / max(abs(eeg_all));

x = eeg_all(1:N);
X_ref = fft(x)/sqrt(N);
F = @(m,a) (1+m).*2.^(-a);

m_vals = linspace(0,1,128);
a_vals = linspace(0,1,128);
[M,A] = meshgrid(m_vals,a_vals);
F_exact = F(M,A);
Y = F_exact(:);

%% ---- Baseline (no interaction)
X_base = [ones(numel(M),1) M(:) A(:)];
coef_base = X_base \ Y;
F_base = reshape(X_base*coef_base,size(M));
err_base = F_base - F_exact;
MRED_base = mean(abs(err_base(:)./F_exact(:)));

%% ---- Proposed_8 (8-bit cross term)
bits = 8;
A_q = round(A*(2^bits-1))/(2^bits-1);
X_prop = [ones(numel(M),1) M(:) A(:) (M(:).*A_q(:))];
coef_prop = X_prop \ Y;
F_prop = reshape(X_prop*coef_prop,size(M));
err_prop = F_prop - F_exact;
MRED_prop = mean(abs(err_prop(:)./F_exact(:)));

stages = log2(N);

sigma_base = MRED_base * sqrt(stages);
sigma_prop = MRED_prop * sqrt(stages);

X_base_fft = X_ref .* (1 + sigma_base*randn(size(X_ref)));
X_prop_fft = X_ref .* (1 + sigma_prop*randn(size(X_ref)));

k = 1:N/2;

figure('Color','w','Position',[100 100 800 600]);

subplot(2,1,1)
plot(k, abs(X_ref(k)),'k','LineWidth',1.6); hold on;
plot(k, abs(X_base_fft(k)),'r--','LineWidth',1);
plot(k, abs(X_prop_fft(k)),'b','LineWidth',1);
grid on; xlim([1 150]);
ylabel('Magnitude');
title('(a) Frequency-Domain Magnitude Comparison');
legend('Exact FFT','Baseline CM','Proposed CM');
set(gca,'FontName','Times New Roman','FontSize',12);

subplot(2,1,2)
plot(k, abs(X_base_fft(k)-X_ref(k)),'r--','LineWidth',1); hold on;
plot(k, abs(X_prop_fft(k)-X_ref(k)),'b','LineWidth',1);
grid on; xlim([1 150]);
xlabel('Frequency Bin');
ylabel('Error Magnitude');
title('(b) Spectral Error Magnitude');
legend('Baseline Error','Proposed Error');
set(gca,'FontName','Times New Roman','FontSize',12);

exportgraphics(gcf,'Frequency_Domain_Comparison.pdf','Resolution',600);