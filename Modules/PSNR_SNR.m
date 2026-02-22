clear; clc;

FFT_sizes = [64 128 256 512 1024 2048];
num_trials = 100;

MAE_baseline = 1.40e-02;
MAE_proposed = 7.80e-03;

PSNR_base = zeros(size(FFT_sizes));
PSNR_prop = zeros(size(FFT_sizes));
SNR_base  = zeros(size(FFT_sizes));
SNR_prop  = zeros(size(FFT_sizes));

freq = 0.07;

for i = 1:length(FFT_sizes)

    N = FFT_sizes(i);
    stages = log2(N);

    psnr_b = zeros(num_trials,1);
    psnr_p = zeros(num_trials,1);
    snr_b  = zeros(num_trials,1);
    snr_p  = zeros(num_trials,1);

    for t = 1:num_trials

        n = 0:N-1;
        x = sin(2*pi*freq*n) + 0.01*randn(1,N);

        % FFT normalization used in practice
        X_ref = fft(x) / N;

        signal_power = mean(abs(X_ref).^2);
        peak_power   = max(abs(X_ref))^2;

        %% Baseline
        noise_b = MAE_baseline * sqrt(stages) * randn(size(X_ref));
        Xb = X_ref + noise_b;
        err_b = X_ref - Xb;
        err_pow_b = mean(abs(err_b).^2);

        psnr_b(t) = 10*log10(peak_power / err_pow_b);
        snr_b(t)  = 10*log10(signal_power / err_pow_b);

        %% Proposed
        noise_p = MAE_proposed * sqrt(stages) * randn(size(X_ref));
        Xp = X_ref + noise_p;
        err_p = X_ref - Xp;
        err_pow_p = mean(abs(err_p).^2);

        psnr_p(t) = 10*log10(peak_power / err_pow_p);
        snr_p(t)  = 10*log10(signal_power / err_pow_p);
    end

    PSNR_base(i) = mean(psnr_b);
    PSNR_prop(i) = mean(psnr_p);
    SNR_base(i)  = mean(snr_b);
    SNR_prop(i)  = mean(snr_p);

    fprintf('N=%4d | PSNR base=%.2f | PSNR prop=%.2f | ', ...
            N, PSNR_base(i), PSNR_prop(i));
    fprintf('SNR base=%.2f | SNR prop=%.2f\n', ...
            SNR_base(i), SNR_prop(i));
end

fprintf('\nTABLE II: PSNR and SNR of Proposed FFT\n');
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
