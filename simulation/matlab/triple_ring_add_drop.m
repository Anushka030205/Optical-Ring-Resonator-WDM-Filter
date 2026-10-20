clc; clear; close all;

%% ==============================
% Wavelength Sweep
%% ==============================
lambda = linspace(1511.2e-9,1511.9e-9,200000);
lambda0 = 1511.54e-9;

c = 3e8;
omega = 2*pi*c./lambda;
omega0 = 2*pi*c/lambda0;

%% ==============================
% Material Parameters
%% ==============================
ng = 4.2;
neff = 2.42;

%% ==============================
% Ring Lengths
%% ==============================
L1 = 100e-6;
L2 = 100e-6;
L3 = 100e-6;

%% ==============================
% Propagation Loss
%% ==============================
loss_dB_cm = 0.5;
alpha = (loss_dB_cm/4.343)*100;

a1_loss = exp(-alpha*L1);

%% ==============================
% Coupling Coefficients
%% ==============================
k_bus  = 0.01;
k_drop = 0.01;
k12    = 0.02;
k23    = 0.02;

tau_scale = c / (ng * L1);

gamma_i  = (1 - a1_loss)*tau_scale;
gamma_e1 = k_bus  * tau_scale;   % Bus coupling
gamma_e2 = k_drop * tau_scale;   % Drop/Add coupling

kappa12  = k12*tau_scale;
kappa23  = k23*tau_scale;

%% ==============================
% ===== ADD-DROP WITH ADD PORT =====
%% ==============================
T_through = zeros(size(lambda));
T_drop    = zeros(size(lambda));
    
s_in  = 1;      % Input from bus
s_add = 0.02;    % Input from ADD port (NEW)

for i = 1:length(lambda)

    delta = omega(i) - omega0;

    % CMT Matrix
    M = [ (1i*delta + gamma_i + gamma_e1 + gamma_e2)   -1i*kappa12    0;
         -1i*kappa12    (1i*delta + gamma_i)       -1i*kappa23;
          0            -1i*kappa23    (1i*delta + gamma_i) ];

    % Input vector (NOW includes ADD PORT)
    B = [ sqrt(2*gamma_e1)*s_in + sqrt(2*gamma_e2)*s_add;
          0;
          0 ];

    % Solve fields
    A = M \ B;
    a1 = A(1);

    % Outputs
    s_through = s_in - sqrt(2*gamma_e1)*a1;
    s_drop    = sqrt(2*gamma_e2)*a1 - s_add;

    % Power
    T_through(i) = abs(s_through)^2;
    T_drop(i)    = abs(s_drop)^2;
end

%% ==============================
% Plot Through & Drop
%% ==============================
figure
plot(lambda*1e9,T_through,'LineWidth',1.5)
hold on
plot(lambda*1e9,T_drop,'LineWidth',1.5)
xlabel('Wavelength (nm)')
ylabel('Transmission')
title('Triple Ring Add-Drop with ADD Port')
legend('Through Port','Drop Port')
grid on

%% ==============================
% ===== Q FACTOR =====
%% ==============================
[peak_val,idx_peak] = max(T_drop);
lambda_res = lambda(idx_peak);

half_max = peak_val/2;

left_idx  = find(T_drop(1:idx_peak) <= half_max,1,'last');
right_idx = find(T_drop(idx_peak:end) <= half_max,1,'first') + idx_peak - 1;

FWHM = lambda(right_idx) - lambda(left_idx);
Q_real = lambda_res / FWHM;

%% ==============================
% FSR
%% ==============================
FSR = lambda_res^2 / (ng * L1);

%% ==============================
% ===== WDM APPLICATION =====
%% ==============================
channel_spacing = 0.4e-9;
num_channels = 5;

channels = lambda_res + (-2:2)*channel_spacing;

input_spectrum = zeros(size(lambda));

for i = 1:length(channels)
    input_spectrum = input_spectrum + ...
        exp(-((lambda - channels(i)).^2)/(2*(0.08e-9)^2));
end

input_spectrum = input_spectrum / max(input_spectrum);

% ADD PORT signal (example: single channel)
add_signal = exp(-((lambda - lambda_res).^2)/(2*(0.08e-9)^2));

% Combined input (bus + add)
output_drop = (input_spectrum + add_signal) .* T_drop;

%% ==============================
% Plot WDM Filtering
%% ==============================
figure
plot(lambda*1e9,input_spectrum,'--','LineWidth',1)
hold on
plot(lambda*1e9,output_drop,'LineWidth',1.5)
xlabel('Wavelength (nm)')
ylabel('Normalized Power')
title('WDM Filtering with ADD Port')
legend('Input WDM','After Filter')
grid on

%% ==============================
% ===== Crosstalk =====
%% ==============================
channel_power = zeros(1,num_channels);

for i = 1:num_channels
    [~,idx] = min(abs(lambda - channels(i)));
    channel_power(i) = output_drop(idx);
end

desired_channel = ceil(num_channels/2);

channel_power(channel_power < 1e-12) = 1e-12;

crosstalk = 10*log10(channel_power / channel_power(desired_channel));

%% ==============================
% ===== PRINT RESULTS =====
%% ==============================
fprintf('\n===== ADD-DROP WITH ADD PORT RESULTS =====\n');
fprintf('Resonance = %.6f nm\n',lambda_res*1e9);
fprintf('FWHM = %.4f pm\n',FWHM*1e12);
fprintf('Q factor = %.3e\n',Q_real);
fprintf('FSR ≈ %.4f nm\n',FSR*1e9);

fprintf('\n===== WDM Crosstalk =====\n');
for i = 1:num_channels
    fprintf('Channel %d Crosstalk = %.2f dB\n',i,crosstalk(i));
end