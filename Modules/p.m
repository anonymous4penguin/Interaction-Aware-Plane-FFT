clear; clc;

FFT_sizes = [64 128 256 512 1024 2048];

MAE_baseline = 1.40e-02;   % baseline CM
MAE_proposed = 7.80e-03;   % interaction-aware CM

alpha = 0.7;   % plane-term contribution
beta  = 0.3;   % cross-term contribution

stage_MSB = [2 4 4 6 6 6 6 8 6 6 8];   % up to 2048 FFT

noise_floor = 1e-12;

%% Load EEG
eeg_all = [];
for k = 0:35
    data = readmatrix(sprintf('s%02d.csv',k));
    eeg_all = [eeg_all; data(:,1)];
end
eeg_all = eeg_all / max(abs(eeg_all));

PSNR_base = zeros(size(FFT_sizes));
PSNR_prop = zeros(size(FFT_sizes));
SNR_base  = zeros(size(FFT_sizes));
SNR_prop  = zeros(size(FFT_sizes));

fprintf('\nEEG-Based PSNR & SNR Evaluation (Correct Model)\n\n');

for i = 1:length(FFT_sizes)

    N = FFT_sizes(i);
    stages = log2(N);

    x = eeg_all(1:N);
    X_ref = fft(x) / sqrt(N);

    signal_power = mean(abs(X_ref).^2);
    peak_power   = max(abs(X_ref))^2;

    %% Baseline
    sigma_base = MAE_baseline * sqrt(stages);
    err_base = sigma_base * randn(size(X_ref));
    P_err_base = max(mean(abs(err_base).^2), noise_floor);

    PSNR_base(i) = 10*log10(peak_power / P_err_base);
    SNR_base(i)  = 10*log10(signal_power / P_err_base);

    %% Proposed (stage-adaptive)
    sigma2 = 0;
    for s = 1:stages
        MSB = stage_MSB(s);
        sigma_stage = MAE_proposed * sqrt( ...
            alpha^2 + beta^2 * 2^(-2*MSB) );
        sigma2 = sigma2 + sigma_stage^2;
    end

    sigma_prop = sqrt(sigma2);
    err_prop = sigma_prop * randn(size(X_ref));
    P_err_prop = max(mean(abs(err_prop).^2), noise_floor);

    PSNR_prop(i) = 10*log10(peak_power / P_err_prop);
    SNR_prop(i)  = 10*log10(signal_power / P_err_prop);

    fprintf('N=%4d | PSNR base=%5.2f | PSNR prop=%5.2f | ', ...
        N, PSNR_base(i), PSNR_prop(i));
    fprintf('SNR base=%5.2f | SNR prop=%5.2f\n', ...
        SNR_base(i), SNR_prop(i));
end
