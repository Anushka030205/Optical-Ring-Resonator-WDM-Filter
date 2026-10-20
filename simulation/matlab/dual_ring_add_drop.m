clc; clear; close all;

%% ======================================================
% Wavelength range
lambda = linspace(1511.2e-9,1511.9e-9,200000);
lambda0 = 1511.54e-9;

c = 3e8;
omega = 2*pi*c./lambda;
omega0 = 2*pi*c/lambda0;

%% ======================================================
% Refractive index
ng   = 4.2;
neff = 2.42;

%% ======================================================
% Ring length
L = 100e-6;

%% ======================================================
% Propagation loss
loss_dB_cm = 2;
alpha = (loss_dB_cm/4.343)*100;

A = exp(-alpha*L);

%% ======================================================
% --------- COUPLING ---------

tau_scale = c / (ng * L);

k_bus  = 0.01;   % input coupling
k_drop = 0.01;   % drop/add coupling
k12    = 0.02;   % inter-ring

gamma_i  = (1 - A)*tau_scale;
gamma_e1 = k_bus  * tau_scale;   % input port
gamma_e2 = k_drop * tau_scale;   % drop/add port
kappa12  = k12    * tau_scale;

%% ======================================================
% --------- ADD-DROP CMT MODEL ---------

T_through = zeros(size(lambda));
T_drop    = zeros(size(lambda));

s_in  = 1;   % input signal
s_add = 0.05;   % add port input (set 0 for now)

for i = 1:length(lambda)

    delta = omega(i) - omega0;

    M = [ (1i*delta + gamma_i + gamma_e1 + gamma_e2)   -1i*kappa12;
         -1i*kappa12                                  (1i*delta + gamma_i) ];

    % BOTH inputs included now
    B = [ sqrt(2*gamma_e1)*s_in + sqrt(2*gamma_e2)*s_add;
          0 ];

    A_vec = M \ B;
    a1 = A_vec(1);

    % OUTPUT PORTS
    s_through = s_in  - sqrt(2*gamma_e1)*a1;
    s_drop    = s_add - sqrt(2*gamma_e2)*a1;

    T_through(i) = abs(s_through)^2;
    T_drop(i)    = abs(s_drop)^2;

end

%% ======================================================
% Plot Through & Drop
figure
plot(lambda*1e9,T_through,'LineWidth',1.5)
hold on
plot(lambda*1e9,T_drop,'LineWidth',1.5)
xlabel('Wavelength (nm)')
ylabel('Transmission')
title('Dual Ring Add-Drop Resonator (With Add Port)')
legend('Through Port','Drop Port')
grid on

%% ======================================================
% ===== Q FACTOR (DROP PORT) =====

[peak_val,idx_peak] = max(T_drop);
lambda_res = lambda(idx_peak);

half_max = peak_val/2;

left_idx  = find(T_drop(1:idx_peak) <= half_max,1,'last');
right_idx = find(T_drop(idx_peak:end) <= half_max,1,'first') + idx_peak - 1;

FWHM = lambda(right_idx) - lambda(left_idx);
Q_real = lambda_res / FWHM;

fprintf('\n========= DUAL RING ADD-DROP (WITH ADD PORT) =========\n');
fprintf('Resonance = %.6f nm\n',lambda_res*1e9);
fprintf('FWHM = %.4f pm\n',FWHM*1e12);
fprintf('Q factor = %.3e\n',Q_real);
%% ======================================================
% ===== FSR CALCULATION =====

FSR = lambda_res^2 / (ng * L);

fprintf('FSR ≈ %.4f nm\n',FSR*1e9);
%% ======================================================
% ================= WDM APPLICATION =====================

channel_spacing = 0.4e-9;
num_channels = 5;

channels = lambda_res + (-2:2)*channel_spacing;

input_spectrum = zeros(size(lambda));

for i = 1:length(channels)
    input_spectrum = input_spectrum + ...
        exp(-((lambda - channels(i)).^2)/(2*(0.08e-9)^2));
end

input_spectrum = input_spectrum / max(input_spectrum);

% Apply ONLY to input port (add port = 0)
output_drop = input_spectrum .* T_drop;

%% ======================================================
% Plot WDM filtering
figure
plot(lambda*1e9,input_spectrum,'--','LineWidth',1)
hold on
plot(lambda*1e9,output_drop,'LineWidth',1.5)
xlabel('Wavelength (nm)')
ylabel('Normalized Power')
title('WDM Filtering using Dual Ring (Add-Drop)')
legend('Input WDM','After Filter')
grid on

%% ======================================================
% ================= CROSSTALK =====================

channel_power = zeros(1,num_channels);

for i = 1:num_channels
    [~,idx] = min(abs(lambda - channels(i)));
    channel_power(i) = output_drop(idx);
end

desired_channel = ceil(num_channels/2);

channel_power(channel_power < 1e-12) = 1e-12;

crosstalk = 10*log10(channel_power / channel_power(desired_channel));

fprintf('\n========= WDM CROSSTALK =========\n');
for i = 1:num_channels
    fprintf('Channel %d Crosstalk = %.2f dB\n',i,crosstalk(i));
end