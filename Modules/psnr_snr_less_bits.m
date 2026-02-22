clear; clc;

FFT_sizes = [256 512 1024];      % use medium sizes for speed
cross_bits = [8 6 5 4 3];        % bits used ONLY in m*a

PSNR = zeros(length(cross_bits), length(FFT_sizes));
SNR  = zeros(length(cross_bits), length(FFT_sizes));

% Load one EEG file (repeatable experiment)
eeg = readmatrix('s00.csv');
x = eeg(:,2)';
x = x - mean(x);                % remove DC

for b = 1:length(cross_bits)
    B = cross_bits(b);

    for i = 1:length(FFT_sizes)
        N = FFT_sizes(i);
        x_seg = x(1:N);

        X_ref = fft(x_seg);

        X_app = approx_fft_cross_bits(x_seg, B);

        err = X_ref - X_app;

        signal_power = mean(abs(X_ref).^2);
        error_power  = mean(abs(err).^2);

        SNR(b,i)  = 10*log10(signal_power / error_power);
        PSNR(b,i) = 10*log10(max(abs(X_ref))^2 / error_power);

        fprintf('Bits=%d | N=%d | PSNR=%.2f | SNR=%.2f\n', ...
                B, N, PSNR(b,i), SNR(b,i));
    end
end

fprintf('\n==== PSNR (dB) ====\n');
disp(array2table(PSNR, ...
    'RowNames', strcat('Bits_', string(cross_bits)), ...
    'VariableNames', strcat('N_', string(FFT_sizes))));

fprintf('\n==== SNR (dB) ====\n');
disp(array2table(SNR, ...
    'RowNames', strcat('Bits_', string(cross_bits)), ...
    'VariableNames', strcat('N_', string(FFT_sizes))));


function X = approx_fft_cross_bits(x, B)

    N = length(x);

    if N == 1
        X = x;
        return;
    end

    Xe = approx_fft_cross_bits(x(1:2:end), B);
    Xo = approx_fft_cross_bits(x(2:2:end), B);

    k = 0:(N/2-1);
    a = k / N;

    % Quantize a to 8 bits (overall angle)
    A8 = round(a * 256) / 256;

    % Reduce bits ONLY for cross-term
    A_cross = floor(A8 * 2^B) / 2^B;

    % Interaction-aware CM (simplified coefficients)
    pm  = 0.95;
    pa  = 0.78;
    pma = 0.32;
    b   = 0.01;

    m = 1;

    Wr = pm*m + pa*(1-A8) + pma*(m.*(1-A_cross)) + b;
    Wi = -(pm*m + pa*A8 + pma*(m.*A_cross) + b);

    W = Wr + 1j*Wi;
    T = W .* Xo;

    X = zeros(1,N);
    X(1:N/2)     = Xe + T;
    X(N/2+1:end) = Xe - T;
end
