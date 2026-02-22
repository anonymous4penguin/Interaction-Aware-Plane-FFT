clear; clc;

FFT_sizes = [64 128 256 512 1024 2048];

% CM MAE values (from YOUR measured experiments)
MAE_baseline = 1.40e-02;   % Plane-fitted CM (paper baseline)
MAE_proposed = 7.80e-03;   % Interaction-aware CM (Stage 4C)

noise_floor = 1e-12;       % Numerical stability

fprintf('Loading EEG data...\n');

eeg_all = [];

for k = 0:35
    fname = sprintf('s%02d.csv', k);
    data  = readmatrix(fname);
    eeg_all = [eeg_all; data(:,1)];
end

fprintf('Total EEG samples loaded: %d\n\n', length(eeg_all));
eeg_all = eeg_all / max(abs(eeg_all));


PSNR_base = zeros(size(FFT_sizes));
PSNR_prop = zeros(size(FFT_sizes));

SNR_base  = zeros(size(FFT_sizes));
SNR_prop  = zeros(size(FFT_sizes));

fprintf('EEG-Based PSNR & SNR Evaluation\n\n');

for i = 1:length(FFT_sizes)

    N = FFT_sizes(i);

    % Take N samples (wrap if needed)
    x = eeg_all(1:N);
    X_ref = fft(x) / sqrt(N);
    stages = log2(N);

    % CM error accumulation model
    noise_scale = sqrt(stages);

    signal_power = mean(abs(X_ref).^2);
    peak_power   = max(abs(X_ref))^2;

    noise_base = MAE_baseline * noise_scale * randn(size(X_ref));
    X_base = X_ref + noise_base;

    err_base = X_ref - X_base;
    error_power_base = mean(abs(err_base).^2);
    error_power_base = max(error_power_base, noise_floor);

    PSNR_base(i) = 10*log10(peak_power / error_power_base);
    SNR_base(i)  = 10*log10(signal_power / error_power_base);

    noise_prop = MAE_proposed * noise_scale * randn(size(X_ref));
    X_prop = X_ref + noise_prop;

    err_prop = X_ref - X_prop;
    error_power_prop = mean(abs(err_prop).^2);
    error_power_prop = max(error_power_prop, noise_floor);

    PSNR_prop(i) = 10*log10(peak_power / error_power_prop);
    SNR_prop(i)  = 10*log10(signal_power / error_power_prop);

    fprintf('N=%4d | PSNR base=%6.2f dB | PSNR prop=%6.2f dB | ', ...
            N, PSNR_base(i), PSNR_prop(i));
    fprintf('SNR base=%6.2f dB | SNR prop=%6.2f dB\n', ...
            SNR_base(i), SNR_prop(i));
end

fprintf('\nTABLE II: PSNR and SNR of Proposed FFT (EEG Signals)\n');

fprintf('FFT length : ');
fprintf('%6d ', FFT_sizes);
fprintf('\n');

fprintf('PSNR base  : ');
fprintf('%6.2f ', PSNR_base);
fprintf('\n');

fprintf('PSNR prop  : ');
fprintf('%6.2f ', PSNR_prop);
fprintf('\n');

fprintf('SNR  base  : ');
fprintf('%6.2f ', SNR_base);
fprintf('\n');

fprintf('SNR  prop  : ');
fprintf('%6.2f ', SNR_prop);
fprintf('\n');
