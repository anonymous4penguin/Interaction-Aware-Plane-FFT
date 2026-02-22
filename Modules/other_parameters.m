clear; clc;

FFT_sizes = [256 512 1024];
bit_list  = [0 2 4 6 8];

eeg_all = [];
for k = 0:35
    data = readmatrix(sprintf('s%02d.csv',k));
    eeg_all = [eeg_all; data(:,1)];
end
eeg_all = eeg_all / max(abs(eeg_all));

F = @(m,a) (1+m).*2.^(-a);
m_vals = linspace(0,1,128);
a_vals = linspace(0,1,128);
[M,A] = meshgrid(m_vals,a_vals);
F_exact = F(M,A);
Y = F_exact(:);

domain_MRED = zeros(size(bit_list));

for b = 1:length(bit_list)

    bits = bit_list(b);

    if bits == 0
        Xreg = [ones(numel(M),1) M(:) A(:)];
        coef = Xreg \ Y;
        F_model = reshape(Xreg*coef,size(M));
    else
        A_q = round(A*(2^bits-1))/(2^bits-1);
        Xreg = [ones(numel(M),1) M(:) A(:) (M(:).*A_q(:))];
        coef = Xreg \ Y;
        F_model = reshape(Xreg*coef,size(M));
    end

    err = F_model - F_exact;
    domain_MRED(b) = mean(abs(err(:)./F_exact(:)));
end

NMSE_all      = zeros(length(bit_list),1);
MaxMagDev_all = zeros(length(bit_list),1);
RSD_all       = zeros(length(bit_list),1);
MaxSpecErr_all= zeros(length(bit_list),1);

N = 1024;
x = eeg_all(1:N);
X_ref = fft(x)/sqrt(N);

signal_energy = sum(abs(X_ref).^2);

stages = log2(N);

for b = 1:length(bit_list)

    % Accumulated multiplicative error variance
    sigma = domain_MRED(b) * sqrt(stages);

    % Inject multiplicative noise
    rng(1); 
    noise = sigma * randn(size(X_ref));

    X_app = X_ref .* (1 + noise);

    err = X_app - X_ref;

    NMSE_all(b) = sum(abs(err).^2) / signal_energy;

    MaxMagDev_all(b) = max(abs(abs(X_app)-abs(X_ref)));

    RSD_all(b) = sqrt(sum((abs(X_app)-abs(X_ref)).^2)) ...
                 / sqrt(sum(abs(X_ref).^2));

    MaxSpecErr_all(b) = max(abs(abs(X_app)-abs(X_ref)) ...
                          ./ (abs(X_ref)+1e-12));
end

fprintf('\n===== Physically-Correct Spectral Metrics =====\n\n');

for b = 1:length(bit_list)

    if bit_list(b)==0
        label = 'Baseline';
    else
        label = sprintf('Proposed_%d',bit_list(b));
    end

    fprintf('%s:\n',label);
    fprintf('NMSE = %.4e\n',NMSE_all(b));
    fprintf('Max Magnitude Deviation = %.4e\n',MaxMagDev_all(b));
    fprintf('Relative Spectral Distortion = %.4e\n',RSD_all(b));
    fprintf('Max Spectral Magnitude Error = %.4e\n\n',MaxSpecErr_all(b));
end