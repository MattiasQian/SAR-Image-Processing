clear all
clc

load sarlab.mat

% Parameters
rho_min    = 10446;
delta_rho  = 4;
theta_deg  = 3;
lambda_c   = 0.0566;
vx         = 131;
fR_over_vx = 2.32;
delta_x    = 0.43;
D          = 1.3;

%% Task 1

%Task 1a
rho_far = rho_min + 349 * delta_rho
L_far   = rho_far * (lambda_c/D)

% Task 1b
M = L_far / delta_x
M = 2* round(M/2)

% Task 1c
rho_37 = rho_min + 36 * delta_rho
L_37   = rho_37 * (lambda_c/D)

%% Task 2

% Task 2a
spectrumplot(signal, 200)

% Task 2b
Baz_theory = (4 * vx / lambda_c) * sind(theta_deg / 2)

%% Task 3

Baz_measured = 165;

% Task 3a
delta_L = lambda_c * rho_far / (2* L_far);

delta_Baz_measured = vx / Baz_measured;

delta_Baz_theoretical = vx / Baz_theory;

delta_theta = lambda_c / (2 * 3 * (pi/180));

delta_D = D / 2;

%%  Task 4

% Task 4a
figure;
plot(abs(oneline).^2);
xlabel('Azimuth sample');
ylabel('Magnitude squared');
title('Squared magnitude of oneline (row 37)');


%% Task 4b                           
m = -M/2 : M/2;                          

%equation (18)
haz = exp(1j * 2*pi * m.^2 / (lambda_c * rho_37 * fR_over_vx^2));

% Apply filter via convolution
oneline_filtered = conv(oneline, haz, 'same');

% Plot
figure;
plot(abs(oneline_filtered).^2);
xlabel('Azimuth sample');
ylabel('Magnitude squared');
title('Filtered oneline (row 37)');

%% Task 4c - Time domain
tic
oneline_td = conv(oneline, haz, 'same');
t_td = toc;

% Frequency domain
tic
N = length(oneline) + length(haz) - 1;
oneline_fd = ifft(fft(oneline, N) .* fft(haz, N));
t_fd = toc;

fprintf('Time domain:      %.4f s\n', t_td);
fprintf('Frequency domain: %.4f s\n', t_fd);
%% Task 5a 
fR_vx_range = 2.31 : 0.01 : 2.34;
resolution = zeros(size(fR_vx_range));

figure;
for k = 1:length(fR_vx_range)
    % Construct filter
    haz_k = exp(1j * 2*pi * m.^2 / (lambda_c * rho_37 * fR_vx_range(k)^2));

    % Apply filter
    oneline_k = conv(oneline, haz_k, 'same');

    % Interpolate 10x
    x_orig = 1:length(oneline_k);
    x_interp = 1:0.1:length(oneline_k);
    oneline_interp = interp1(x_orig, abs(oneline_k), x_interp);
    mag_sq = oneline_interp.^2;

    % Find peak and 3dB width
    [peak, idx] = max(mag_sq);
    half_power = peak / 2;

    left  = find(mag_sq(1:idx) < half_power, 1, 'last');
    right = find(mag_sq(idx:end) < half_power, 1, 'first') + idx - 1;

    resolution(k) = (right - left) * delta_x / 10;

    % Close-up plot of the reflector
    subplot(2, 2, k);
    window = (idx-200):(idx+200);
    plot(x_interp(window) * delta_x, mag_sq(window));
    xlabel('Azimuth position (m)');
    ylabel('Magnitude squared');
    title(sprintf('f_R/v_x = %.2f m^{-1}, \\delta = %.2f m', ...
        fR_vx_range(k), resolution(k)));
    grid on;
end

% Print table of results
fprintf('\n  fR/vx [m^-1]   |  3-dB width [m]\n');
for k = 1:length(fR_vx_range)
    fprintf('     %.2f       |     %.2f\n', fR_vx_range(k), resolution(k));
end

[min_res, idx_opt] = min(resolution);
fR_vx_opt = fR_vx_range(idx_opt);
fprintf('\nOptimal fR/vx = %.2f m^-1, min 3-dB width = %.2f m\n', fR_vx_opt, min_res);

%% Task 6

% Task 6a
figure
imageplot(signal);
title('Compressed data');

%% Task 6b

% Build filter matrix
Hfilter = zeros(350, M+1);
m = -M/2 : M/2;

for row = 1:350
    rho_row = rho_min + (row-1) * delta_rho;
    Hfilter(row, :) = exp(1j * 2*pi * m.^2 / (lambda_c * rho_row * fR_vx_opt^2));
end

image = sarfilter(Hfilter, signal);

% Display image
figure;
imageplot(image);
title('Filtered SAR image');

%% Task 7 

% Task 7a
spectrumplot(image, 200);  % average over same number of rows as Task 2a
title('Azimuth power spectrum of focused image');

%% Task 7b
filterplot(Hfilter(1, :))

%% Task 8

%% Task 8a
w = ones(1, 9) / 9;                          % 9-look filter
image_multi = sqrt(conv2(abs(image).^2, w, 'same'));

% Display the 9-look image
figure;
imageplot(image_multi);
title('9-look SAR image');

%% Task 8 b&c
oneline_multi = image_multi(37, :);     % row 37 from multi-look image

% Interpolate 10x
x_orig = 1:length(oneline_multi);
x_interp = 1:0.1:length(oneline_multi);
oneline_interp = interp1(x_orig, oneline_multi, x_interp);

% Find peak around column 2180 and measure 3-dB width
[peak, idx] = max(oneline_interp);
half_power = peak / 2;

left  = find(oneline_interp(1:idx) < half_power, 1, 'last');
right = find(oneline_interp(idx:end) < half_power, 1, 'first') + idx - 1;

resolution_9l = (right - left) * delta_x / 10;
fprintf('Measured 9-look resolution: %.2f m\n', resolution_9l);
fprintf('Theoretical 9-look resolution: %.2f m\n', 9 * delta_D);