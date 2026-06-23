%% Startup cleanup
clearvars; close all;

%% Material properties
% Properties of the CFRP T800/913
rho = 1550.0;  % [kg/m^3]
C_CFRP = [154.0, 3.7,    3.7,    0,      0,      0;  % [GPa]
          3.7,   9.5,    5.2,    0,      0,      0;
          3.7,   5.2,    9.5,    0,      0,      0;
          0,     0,      0,      2.15,   0,      0;
          0,     0,      0,      0,      4.2,    0;
          0,     0,      0,      0,      0,      4.2] * 1e9;

% Properties of the aluminium substrate
rho_sub = 2700.0; % [kg/m^3]
C_sub = [107.5   54.59     54.59     0       0       0; % [GPa]
         54.59   107.5     54.59     0       0       0;
         54.59   54.59     107.5     0       0       0;
         0       0         0         26.45   0       0;
         0       0         0         0       26.45   0;
         0       0         0         0       0       26.45] * 1e9;

lambda_sub = C_sub(1,2);
mu_sub = C_sub(5,5);

ct_sub = sqrt(mu_sub/rho_sub);
cl_sub = sqrt((lambda_sub + 2*mu_sub)/rho_sub);

%% Root filter parameters
k_abs_max = 1e7;            % Maximum absolute wavenumber value
att_min_tol = -1e-6;        % With exp(i*k*x-i*w*t), Im(k)>=0 is attenuating
att_max = 15e3;             % Maximum attenuation value
qep_res_tol = 1e-6;         % QEP residual tolerance
interface_res_tol = 1e-5;   % Analytical substrate-interface residual tolerance
kh_min = 0.116;             % Minimum kh for fully leaky branch

%% Multilayer properties and general parameters
% Each array/cell value corresponds to one layer from top to bottom.
h_layer = 0.5e-3; % [m]
densities = [rho rho rho rho];
Cs = {C_CFRP C_CFRP C_CFRP C_CFRP};
orientations = [0 0 0 0];       % [degrees]
thicknesses = [h_layer h_layer h_layer h_layer];
Ns = [70 70 70 70];             % Collocation points for each plate layer

num_of_layers = length(thicknesses);
total_h = sum(thicknesses);
propagation_angle = 0.0;        % [degrees]

N_sub = 70;                     % Substrate collocation points
freq_limit = 2e6;               % [Hz]
PV_limit = 12000.0;             % [m/s]
WN_limit = 8000.0;              % [1/m]
F_amount = 150;


%% Transform stiffness matrices relative to propagation direction
beta = orientations - propagation_angle;
c = cell(1, num_of_layers);
for n = 1:num_of_layers
    c{n} = transform_stiffness_matrix(Cs{n}, beta(n));
end

%% Half-space map parameters
% Real zeta gives evanescent rational-Chebyshev mapping.
% Imaginary zeta gives complex continuation for leaky radiation.
% zeta_phi controls the longitudinal potential; zeta_psi controls the shear potential.
multiplier = 20;
cases(1).name = 'evanescent';
cases(1).zeta_phi = multiplier * total_h;
cases(1).zeta_psi = multiplier * total_h;

cases(2).name = 'shear-leaky';
cases(2).zeta_phi = multiplier * total_h;
cases(2).zeta_psi = 1i * multiplier * total_h;

cases(3).name = 'fully-leaky';
cases(3).zeta_phi = 1i * multiplier * total_h;
cases(3).zeta_psi = 1i * multiplier * total_h;

%% Multilayer plate setup
DM1 = cell(1, num_of_layers);
DM2 = cell(1, num_of_layers);
z = cell(1, num_of_layers);
z_mm = cell(1, num_of_layers);
I = cell(1, num_of_layers);

for n = 1:num_of_layers
    [s, DM] = chebdif(Ns(n), 2);

    z_min = total_h/2 - sum(thicknesses(1:n));
    z_max = total_h/2 - sum(thicknesses(1:(n-1)));
    z_c = 0.5 * (z_min + z_max);

    scale = 2.0 / (z_max - z_min);
    DM1{n} = DM(:,:,1) * scale;
    DM2{n} = DM(:,:,2) * scale^2;

    z{n} = 0.5 * (z_max - z_min) * s + z_c;
    z_mm{n} = z{n} * 1e3;
    I{n} = eye(Ns(n));
end

L0_layers = cell(1, num_of_layers);
L1_layers = cell(1, num_of_layers);
L2_layers = cell(1, num_of_layers);
S0_layers = cell(1, num_of_layers);
S1_layers = cell(1, num_of_layers);
M_layers = cell(1, num_of_layers);

for n = 1:num_of_layers
    [L0_layers{n}, L1_layers{n}, L2_layers{n}] = get_Lamb_L_coeffs(c{n}, I{n}, DM1{n}, DM2{n});
    [S0_layers{n}, S1_layers{n}] = get_Lamb_S_coeffs(c{n}, I{n}, DM1{n});
    M_layers{n} = get_Lamb_M_matrix(Ns(n), densities(n));
end

%% Frequency sweep setup
freqs = linspace(1.0, freq_limit, F_amount);
omegas = 2.0 * pi * freqs;

%% Substrate cases and calculation
freq_all = [];
cph_all = [];
att_all = [];
k_all = [];
res_all = [];
bc_all = [];
case_all = [];
U_modes = [];

p_size = 2 * sum(Ns);
z_interface = -total_h/2;

t0 = tic;
for case_i = 1:numel(cases)
    zeta_phi_sub = cases(case_i).zeta_phi;
    zeta_psi_sub = cases(case_i).zeta_psi;

    [~, ~, DM1_phi_sub, DM2_phi_sub] = chebdif_sub(N_sub, z_interface, zeta_phi_sub);
    [~, ~, DM1_psi_sub, DM2_psi_sub] = chebdif_sub(N_sub, z_interface, zeta_psi_sub);
    I_sub = eye(N_sub);

    M_sub = get_Lamb_M_matrix(N_sub, rho_sub);
    [L0_sub, L1_sub, L2_sub] = get_Lamb_potential_L_coeffs(lambda_sub, mu_sub, I_sub, DM2_phi_sub, DM2_psi_sub);
    [U0_sub, U1_sub, U2_sub] = get_Lamb_potential_U_coeffs(I_sub, DM1_phi_sub, DM1_psi_sub);
    [T0_sub, T1_sub, T2_sub] = get_Lamb_potential_T_coeffs(lambda_sub, mu_sub, I_sub, DM1_phi_sub, DM2_phi_sub, DM1_psi_sub, DM2_psi_sub);

    [L0_global, L1_global, L2_global, M_global] = assemble_Lamb_ML_substrate_QEP_coeffs( ...
        L0_layers, L1_layers, L2_layers, S0_layers, S1_layers, I, Ns, M_layers, ...
        L0_sub, L1_sub, L2_sub, M_sub, U0_sub, U1_sub, U2_sub, T0_sub, T1_sub, T2_sub, N_sub);

    freq_plot = [];
    cph_plot = [];
    att_plot = [];
    k_plot = [];
    res_plot = [];
    bc_res_plot = [];

    fprintf('Starting multilayer Lamb QEP calculation, substrate case: %s\n', cases(case_i).name);

    for i = 1:numel(omegas)
        omega = omegas(i);

        if mod(i-1, 25) == 0
            pct = 100.0 * (i-1) / numel(omegas);
            fprintf('  %.1f %%\n', pct);
        end

        A0 = L0_global - (omega^2) * M_global;
        A1 = L1_global;
        A2 = L2_global;

        [kvals, U] = solve_QEP_quadeig(A0, A1, A2);

        for m = 1:numel(kvals)
            k = kvals(m);

            % Keep only waves propagating in the positive x direction.
            if real(k) <= 0
                continue;
            end

            if ~(isfinite(real(k)) && isfinite(imag(k)))
                continue;
            end

            alpha = imag(k);
            if alpha < att_min_tol
                continue;
            end
            if alpha > att_max
                continue;
            end
            if abs(k) > k_abs_max
                continue;
            end

            cph = omega / real(k);
            if ~(isfinite(cph) && cph > 0)
                continue;
            end

            % Select the proper substrate continuation according to phase velocity.
            if case_i == 1
                if cph >= ct_sub
                    continue;
                end
            elseif case_i == 2
                if cph < ct_sub || cph >= cl_sub
                    continue;
                end
            elseif case_i == 3
                if cph < cl_sub
                    continue;
                end
                kh_val = real(k) * total_h;
                if ~(isfinite(kh_val) && kh_val >= kh_min)
                    continue;
                end
            end

            u = U(:, m);
            if any(~isfinite(u))
                continue;
            end

            res = qep_residual(A0, A1, A2, k, u);
            if res > qep_res_tol
                continue;
            end

            bc_res = interface_residual_ML( ...
                k, omega, u, Ns, N_sub, p_size, ...
                lambda_sub, mu_sub, rho_sub, zeta_phi_sub, zeta_psi_sub, ...
                S0_layers{end}, S1_layers{end});

            if ~(isfinite(bc_res) && bc_res <= interface_res_tol)
                continue;
            end

            freq_plot(end+1,1) = omega / (2*pi); %#ok<AGROW>
            cph_plot(end+1,1) = real(cph); %#ok<AGROW>
            att_plot(end+1,1) = imag(k); %#ok<AGROW>
            k_plot(end+1,1) = k; %#ok<AGROW>
            res_plot(end+1,1) = res; %#ok<AGROW>
            bc_res_plot(end+1,1) = bc_res; %#ok<AGROW>
            U_modes(:,end+1) = u; %#ok<AGROW>

            % Optional displacement-profile check for accepted roots.
            % [Ux_layers, Uz_layers] = get_Lamb_ML_displacement_profile(u(1:p_size), Ns);
            % figure(100); clf; hold on;
            % for layer_i = 1:num_of_layers
            %     Ux_plot = real(Ux_layers{layer_i});
            %     Uz_plot = real(Uz_layers{layer_i});
            %     if max(abs(Ux_plot)) > 0
            %         Ux_plot = Ux_plot / max(abs(Ux_plot));
            %     end
            %     if max(abs(Uz_plot)) > 0
            %         Uz_plot = Uz_plot / max(abs(Uz_plot));
            %     end
            %     plot(Ux_plot, z_mm{layer_i}, '-o', 'DisplayName', sprintf('Re(Ux), layer %d', layer_i));
            %     plot(Uz_plot, z_mm{layer_i}, '-^', 'DisplayName', sprintf('Re(Uz), layer %d', layer_i));
            % end
            % xlabel('Normalized amplitude'); ylabel('z (mm)');
            % title(sprintf('f = %.2f kHz, c_p = %.1f m/s', omega/(2*pi)/1e3, cph));
            % legend('Location','best'); grid on; drawnow;
        end
    end

    freq_all = [freq_all; freq_plot]; %#ok<AGROW>
    cph_all = [cph_all; cph_plot]; %#ok<AGROW>
    att_all = [att_all; att_plot]; %#ok<AGROW>
    k_all = [k_all; k_plot]; %#ok<AGROW>
    res_all = [res_all; res_plot]; %#ok<AGROW>
    bc_all = [bc_all; bc_res_plot]; %#ok<AGROW>
    case_all = [case_all; case_i * ones(size(freq_plot))]; %#ok<AGROW>

    fprintf('Case %s accepted roots: %d\n', cases(case_i).name, numel(k_plot));
end

fprintf('  100.0 %%\n');
fprintf('Multilayer QEP calculation took %.3f s\n', toc(t0));
fprintf('Accepted roots: %d\n', numel(k_all));

%% Plots
k_real = real(k_all);

freq_line_kHz = linspace(0, freq_limit*1e-3, 500).';
omega_line = 2*pi*freq_line_kHz*1e3;
k_CT_Al = omega_line / ct_sub;
k_CL_Al = omega_line / cl_sub;

% Wavenumber versus frequency plot
figure('Color','w'); hold on;
idx = case_all == 1;
plot(freq_all(idx)*1e-3, k_real(idx), '.r', 'MarkerSize', 4);
idx = case_all == 2;
plot(freq_all(idx)*1e-3, k_real(idx), '.b', 'MarkerSize', 4);
idx = case_all == 3;
plot(freq_all(idx)*1e-3, k_real(idx), '.k', 'MarkerSize', 4);
plot(freq_line_kHz, k_CT_Al, '--k', 'LineWidth', 1.0);
plot(freq_line_kHz, k_CL_Al, ':k',  'LineWidth', 1.0);
xlim([0, freq_limit*1e-3]);
ylim([0, WN_limit]);
xlabel('Frequency (kHz)');
ylabel('Wavenumber Re(k_x) (1/m)');
legend('evanescent','shear-leaky','fully-leaky', ...
       'c_T Al boundary','c_L Al boundary', ...
       'Location','best');
grid on;

% Phase velocity versus frequency plot
figure('Color','w'); hold on;
idx = case_all == 1;
plot(freq_all(idx)*1e-3, cph_all(idx), '.r', 'MarkerSize', 4);
idx = case_all == 2;
plot(freq_all(idx)*1e-3, cph_all(idx), '.b', 'MarkerSize', 4);
idx = case_all == 3;
plot(freq_all(idx)*1e-3, cph_all(idx), '.k', 'MarkerSize', 4);
yline(ct_sub, '--k', 'c_T Al');
yline(cl_sub, '--k', 'c_L Al');
xlim([0, freq_limit*1e-3]);
ylim([0, PV_limit]);
xlabel('Frequency (kHz)');
ylabel('Phase velocity (m/s)');
legend('evanescent','shear-leaky','fully-leaky','Location','best');
grid on;

% All accepted roots, without substrate-case coloring
figure('Color','w');
plot(freq_all*1e-3, cph_all, '.k', 'MarkerSize', 4);
yline(ct_sub, '--k', 'c_T Al');
yline(cl_sub, '--k', 'c_L Al');
xlim([0, freq_limit*1e-3]);
ylim([0, PV_limit]);
xlabel('Frequency (kHz)');
ylabel('Phase velocity (m/s)');
grid on;

% Attenuation versus frequency plot
figure('Color','w'); hold on;
idx = case_all == 1;
plot(freq_all(idx)*1e-3, att_all(idx), '.r', 'MarkerSize', 4, 'DisplayName', 'evanescent');
idx = case_all == 2;
plot(freq_all(idx)*1e-3, att_all(idx), '.b', 'MarkerSize', 4, 'DisplayName', 'shear-leaky');
idx = case_all == 3;
plot(freq_all(idx)*1e-3, att_all(idx), '.k', 'MarkerSize', 4, 'DisplayName', 'fully-leaky');
xlim([0, freq_limit*1e-3]);
ylim([0, att_max]);
xlabel('Frequency (kHz)');
ylabel('Attenuation Im(k_x) (1/m)');
legend('Location','best');
grid on; box on;

%% Functions
function [kvals, U] = solve_QEP_quadeig(A0, A1, A2)
% Balanced quadeig solver for:
%     (A2*k^2 + A1*k + A0) U = 0
% with k = k0*lambda.

n = size(A0, 1);
k0 = 1e4;

B0 = complex(A0);
B1 = complex(A1) * k0;
B2 = complex(A2) * k0^2;

row_scale = max([abs(B0), abs(B1), abs(B2)], [], 2);
zero_rows = find(row_scale == 0);
if ~isempty(zero_rows)
    error('Zero rows found in QEP coefficients. First zero row: %d', zero_rows(1));
end

R = spdiags(1 ./ row_scale, 0, n, n);
B0 = R * B0;
B1 = R * B1;
B2 = R * B2;

try
    [U, D] = quadeig(B2, B1, B0);
    if isvector(D)
        lambda_vals = D(:);
    else
        lambda_vals = diag(D);
    end
    kvals = k0 * lambda_vals;
catch ME
    fprintf('\nBalanced quadeig failed with this error:\n');
    fprintf('%s\n', ME.message);
    rethrow(ME);
end
end

function [x, DM_out] = chebdif(N, M)
I_mat = eye(N);
L = logical(I_mat);

n1 = floor(N/2);
n2 = floor((N+1)/2);

k = (0:N-1)';
th = k * pi / (N - 1);

x = sin(pi * ((N-1):-2:-(N-1))' / (2 * (N - 1)));

T = repmat(th/2, 1, N);
DX = 2.0 * sin(T' + T) .* sin(T' - T);
DX = [DX(1:n1,:); -flipud(fliplr(DX(1:n2,:)))];
DX(L) = 1.0;

c_vec = (-1.0).^k;
C = toeplitz(c_vec);
C(1,:) = 2*C(1,:);
C(end,:) = 2*C(end,:);
C(:,1) = C(:,1)/2;
C(:,end) = C(:,end)/2;

Z = 1.0 ./ DX;
Z(L) = 0.0;

D_mat = eye(N);
DM_out = zeros(N, N, M);

for ell = 1:M
    D_mat = ell * Z .* (C .* repmat(diag(D_mat), 1, N) - D_mat);
    D_mat(1:N+1:end) = -sum(D_mat, 2);
    DM_out(:,:,ell) = D_mat;
end
end

function [s_sub, z_sub, D1_sub, D2_sub] = chebdif_sub(Nsub, z_interface, zeta)
% Chebyshev differentiation matrices for a bottom semi-infinite substrate.
% Physical domain:
%       z in (-inf, z_interface]
% Mapping:
%       z(s) = z_interface - zeta * (1 - s)/(1 + s)
% With this chebdif ordering:
%       s(1)    = +1  -> z(1)    = z_interface
%       s(Nsub) = -1  -> z(Nsub) = -inf

[s_sub, DM] = chebdif(Nsub, 2);
Ds1 = DM(:,:,1);
Ds2 = DM(:,:,2);

den = 1.0 + s_sub;
finite_nodes = abs(den) > 100*eps;
z_sub = complex(nan(Nsub,1));
z_sub(finite_nodes) = z_interface ...
    - zeta * (1.0 - s_sub(finite_nodes)) ./ (1.0 + s_sub(finite_nodes));
z_sub(~finite_nodes) = NaN;

% Chain-rule factors: d/dz = (ds/dz) d/ds
dsdz   = (1.0 + s_sub).^2 / (2.0*zeta);
d2sdz2 = (1.0 + s_sub).^3 / (2.0*zeta^2);

dsdz(~finite_nodes) = 0.0;
d2sdz2(~finite_nodes) = 0.0;

D1_sub = diag(dsdz) * Ds1;
D2_sub = diag(dsdz.^2) * Ds2 + diag(d2sdz2) * Ds1;
end

function c = transform_stiffness_matrix(C, beta)
% Rotation about the thickness-normal convention used in your previous code.
c = zeros(6, 6);
s = sind(beta);
g = cosd(beta);

c(1, 1) = C(1,1)*g^4 + C(2,2)*s^4 + 2*(C(1,2)+2*C(6,6))*s^2*g^2;
c(1, 2) = (C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2 + C(1,2);
c(1, 3) = C(1,3)*g^2 + C(2,3)*s^2;
c(1, 6) = (C(1,2)+2*C(6,6)-C(1,1))*s*g^3 + (C(2,2)-C(1,2)-2*C(6,6))*g*s^3;
c(2, 2) = C(1,1)*s^4 + C(2,2)*g^4 + 2*(C(1,2)+2*C(6,6))*s^2*g^2;
c(2, 3) = C(2,3)*g^2 + C(1,3)*s^2;
c(2, 6) = (C(1,2)+2*C(6,6)-C(1,1))*g*s^3 + (C(2,2)-C(1,2)-2*C(6,6))*s*g^3;
c(3, 3) = C(3,3);
c(3, 6) = (C(2,3)-C(1,3))*s*g;
c(4, 4) = C(4,4)*g^2 + C(5,5)*s^2;
c(4, 5) = (C(4,4)-C(5,5))*s*g;
c(5, 5) = C(5,5)*g^2 + C(4,4)*s^2;
c(6, 6) = C(6,6) + (C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2;

c(2,1)=c(1,2); c(3,1)=c(1,3); c(6,1)=c(1,6);
c(3,2)=c(2,3); c(6,2)=c(2,6);
c(6,3)=c(3,6); c(5,4)=c(4,5);
end

function [L0, L1, L2] = get_Lamb_L_coeffs(c, I, D, D2)
% Displacement order: [Ux; Uz]
%
% Polynomial coefficients for:
%     L(k) = L0 + k*L1 + k^2*L2

L11_0 = c(5,5) * D2;
L13_0 = c(3,5) * D2;
L31_0 = L13_0;
L33_0 = c(3,3) * D2;

L0 = [L11_0, L13_0;
      L31_0, L33_0];

L11_1 = 2i * c(1,5) * D;
L13_1 = 1i * (c(1,3) + c(5,5)) * D;
L31_1 = L13_1;
L33_1 = 2i * c(3,5) * D;

L1 = [L11_1, L13_1;
      L31_1, L33_1];

L11_2 = -c(1,1) * I;
L13_2 = -c(1,5) * I;
L31_2 = L13_2;
L33_2 = -c(5,5) * I;

L2 = [L11_2, L13_2;
      L31_2, L33_2];
end

function [S0, S1] = get_Lamb_S_coeffs(c, I, D)
% Lamb traction operator on a horizontal z = const boundary.
% Displacement order: [Ux; Uz]
%
% Row block 1: sigma_zz
% Row block 2: sigma_xz

Szz_Ux_0 = c(3,5) * D;
Szz_Uz_0 = c(3,3) * D;
Sxz_Ux_0 = c(5,5) * D;
Sxz_Uz_0 = c(3,5) * D;

S0 = [Szz_Ux_0, Szz_Uz_0;
      Sxz_Ux_0, Sxz_Uz_0];

Szz_Ux_1 = 1i * c(1,3) * I;
Szz_Uz_1 = 1i * c(3,5) * I;
Sxz_Ux_1 = 1i * c(1,5) * I;
Sxz_Uz_1 = 1i * c(5,5) * I;

S1 = [Szz_Ux_1, Szz_Uz_1;
      Sxz_Ux_1, Sxz_Uz_1];
end

function [L0, L1, L2] = get_Lamb_potential_L_coeffs(lambda_sub, mu_sub, I, D2_phi, D2_psi)
% Aluminium substrate potentials order: [phi; psi]
Z = zeros(size(I));

L0 = [(lambda_sub + 2.0*mu_sub) * D2_phi, Z;
      Z,                                  mu_sub * D2_psi];

L1 = zeros(2*size(I,1), 2*size(I,1));

L2 = [-(lambda_sub + 2.0*mu_sub) * I, Z;
      Z,                                -mu_sub * I];
end

function [U0, U1, U2] = get_Lamb_potential_U_coeffs(I, D_phi, D_psi)
% Aluminium substrate displacement reconstruction from potentials.
% Potential order: [phi; psi]
%
% ux = i*k*phi - D(psi)
% uz = D(phi) + i*k*psi

Z = zeros(size(I));

Ux_0 = [Z, -D_psi];
Uz_0 = [D_phi, Z];
U0 = [Ux_0;
      Uz_0];

Ux_1 = [1i*I, Z];
Uz_1 = [Z, 1i*I];
U1 = [Ux_1;
      Uz_1];

U2 = zeros(2*size(I,1), 2*size(I,1));
end

function [T0, T1, T2] = get_Lamb_potential_T_coeffs(lambda_sub, mu_sub, I, D_phi, D2_phi, D_psi, D2_psi)
% Aluminium substrate tractions from potentials.
% Potential order: [phi; psi]

Z = zeros(size(I));

Szz_0 = [(lambda_sub + 2.0*mu_sub) * D2_phi, Z];
Sxz_0 = [Z, -mu_sub * D2_psi];
T0 = [Szz_0;
      Sxz_0];

Szz_1 = [Z, 2i * mu_sub * D_psi];
Sxz_1 = [2i * mu_sub * D_phi, Z];
T1 = [Szz_1;
      Sxz_1];

Szz_2 = [-lambda_sub * I, Z];
Sxz_2 = [Z, -mu_sub * I];
T2 = [Szz_2;
      Sxz_2];
end

function M = get_Lamb_M_matrix(N, density)
M = -density * eye(N * 2);
end

function [A0, A1, A2, M] = assemble_Lamb_ML_substrate_QEP_coeffs( ...
    L0_layers, L1_layers, L2_layers, S0_layers, S1_layers, Is, Ns, M_layers, ...
    L0_sub, L1_sub, L2_sub, M_sub, U0_sub, U1_sub, U2_sub, T0_sub, T1_sub, T2_sub, N_sub)
% Assemble a decoupled Lamb multilayer plate bonded to an isotropic
% half-space represented by scalar/vector potentials.
%
% Plate layer DOF order inside each layer:
%     [Ux_1 ... Ux_N, Uz_1 ... Uz_N]
% Layers are stacked from top to bottom. Substrate potential order is:
%     [phi_1 ... phi_Nsub, psi_1 ... psi_Nsub]

num_of_layers = length(Ns);
p_size = 2 * sum(Ns);
s_size = 2 * N_sub;
total_size = p_size + s_size;

A0 = complex(zeros(total_size, total_size));
A1 = complex(zeros(total_size, total_size));
A2 = complex(zeros(total_size, total_size));
M  = complex(zeros(total_size, total_size));

% Bulk layer equations.
for n = 1:num_of_layers
    cols = layer_cols(n, Ns);
    A0(cols, cols) = complex(L0_layers{n});
    A1(cols, cols) = complex(L1_layers{n});
    A2(cols, cols) = complex(L2_layers{n});
    M(cols, cols)  = complex(M_layers{n});
end

% Bulk substrate potential equations.
sub_cols = p_size + (1:s_size);
A0(sub_cols, sub_cols) = complex(L0_sub);
A1(sub_cols, sub_cols) = complex(L1_sub);
A2(sub_cols, sub_cols) = complex(L2_sub);
M(sub_cols, sub_cols)  = complex(M_sub);

% 1. Top free surface of the first layer.
N1 = Ns(1);
cols_1 = layer_cols(1, Ns);

row = layer_offset(1, Ns) + 1;
[A0, A1, A2, M] = clear_QEP_row(A0, A1, A2, M, row);
A0(row, cols_1) = S0_layers{1}(1, :);
A1(row, cols_1) = S1_layers{1}(1, :);

row = layer_offset(1, Ns) + N1;
[A0, A1, A2, M] = clear_QEP_row(A0, A1, A2, M, row);
A0(row, cols_1) = S0_layers{1}(N1 + 1, :);
A1(row, cols_1) = S1_layers{1}(N1 + 1, :);

% 2. Internal multilayer interfaces.
for n = 1:(num_of_layers - 1)
    N_left = Ns(n);
    N_right = Ns(n + 1);
    cols_left = layer_cols(n, Ns);
    cols_right = layer_cols(n + 1, Ns);

    % Traction continuity at bottom of layer n / top of layer n+1.
    row = layer_offset(n, Ns) + N_left + 1;
    [A0, A1, A2, M] = clear_QEP_row(A0, A1, A2, M, row);
    A0(row, cols_left)  = S0_layers{n}(N_left, :);
    A1(row, cols_left)  = S1_layers{n}(N_left, :);
    A0(row, cols_right) = -S0_layers{n+1}(1, :);
    A1(row, cols_right) = -S1_layers{n+1}(1, :);

    row = layer_offset(n, Ns) + 2*N_left;
    [A0, A1, A2, M] = clear_QEP_row(A0, A1, A2, M, row);
    A0(row, cols_left)  = S0_layers{n}(2*N_left, :);
    A1(row, cols_left)  = S1_layers{n}(2*N_left, :);
    A0(row, cols_right) = -S0_layers{n+1}(N_right + 1, :);
    A1(row, cols_right) = -S1_layers{n+1}(N_right + 1, :);

    % Displacement continuity at the same interface.
    row = layer_offset(n + 1, Ns) + 1;
    [A0, A1, A2, M] = clear_QEP_row(A0, A1, A2, M, row);
    A0(row, cols_left)  = get_Lamb_interface_BC(Is{n}(N_left, :), N_left, 1);
    A0(row, cols_right) = get_Lamb_interface_BC(-Is{n+1}(1, :), N_right, 1);

    row = layer_offset(n + 1, Ns) + N_right;
    [A0, A1, A2, M] = clear_QEP_row(A0, A1, A2, M, row);
    A0(row, cols_left)  = get_Lamb_interface_BC(Is{n}(N_left, :), N_left, 2);
    A0(row, cols_right) = get_Lamb_interface_BC(-Is{n+1}(1, :), N_right, 2);
end

% 3. Bottom interface between the last plate layer and the substrate.
n = num_of_layers;
NL = Ns(n);
cols_last = layer_cols(n, Ns);
s_i_top = [1, N_sub + 1];
s_i_bot = [N_sub, 2*N_sub];

% Traction continuity.
row = layer_offset(n, Ns) + NL + 1;
[A0, A1, A2, M] = clear_QEP_row(A0, A1, A2, M, row);
A0(row, cols_last) = S0_layers{n}(NL, :);
A1(row, cols_last) = S1_layers{n}(NL, :);
A0(row, sub_cols) = -T0_sub(s_i_top(1), :);
A1(row, sub_cols) = -T1_sub(s_i_top(1), :);
A2(row, sub_cols) = -T2_sub(s_i_top(1), :);

row = layer_offset(n, Ns) + 2*NL;
[A0, A1, A2, M] = clear_QEP_row(A0, A1, A2, M, row);
A0(row, cols_last) = S0_layers{n}(2*NL, :);
A1(row, cols_last) = S1_layers{n}(2*NL, :);
A0(row, sub_cols) = -T0_sub(s_i_top(2), :);
A1(row, sub_cols) = -T1_sub(s_i_top(2), :);
A2(row, sub_cols) = -T2_sub(s_i_top(2), :);

% Displacement continuity.
row = p_size + s_i_top(1);
[A0, A1, A2, M] = clear_QEP_row(A0, A1, A2, M, row);
A0(row, cols_last) = get_Lamb_interface_BC(Is{n}(NL, :), NL, 1);
A0(row, sub_cols) = -U0_sub(s_i_top(1), :);
A1(row, sub_cols) = -U1_sub(s_i_top(1), :);
A2(row, sub_cols) = -U2_sub(s_i_top(1), :);

row = p_size + s_i_top(2);
[A0, A1, A2, M] = clear_QEP_row(A0, A1, A2, M, row);
A0(row, cols_last) = get_Lamb_interface_BC(Is{n}(NL, :), NL, 2);
A0(row, sub_cols) = -U0_sub(s_i_top(2), :);
A1(row, sub_cols) = -U1_sub(s_i_top(2), :);
A2(row, sub_cols) = -U2_sub(s_i_top(2), :);

% 4. Substrate condition at infinity.
for j = 1:2
    row = p_size + s_i_bot(j);
    [A0, A1, A2, M] = clear_QEP_row(A0, A1, A2, M, row);
    A0(row, row) = 1.0;
end
end

function offset = layer_offset(n, Ns)
offset = 2 * sum(Ns(1:(n-1)));
end

function cols = layer_cols(n, Ns)
cols = (layer_offset(n, Ns) + 1):(layer_offset(n, Ns) + 2*Ns(n));
end

function [A0, A1, A2, M] = clear_QEP_row(A0, A1, A2, M, row)
A0(row, :) = 0.0;
A1(row, :) = 0.0;
A2(row, :) = 0.0;
M(row, :) = 0.0;
end

function I_bc = get_Lamb_interface_BC(I_n, N, Pos)
I_bc = zeros(1, 2*N);
I_bc(1, (N*(Pos-1) + 1):N*Pos) = I_n;
end

function r = qep_residual(A0, A1, A2, k, u)
num = norm((A2*k^2 + A1*k + A0) * u);
den = (norm(A2)*abs(k)^2 + norm(A1)*abs(k) + norm(A0)) * norm(u);
if den == 0
    r = inf;
else
    r = num / den;
end
end

function bc_res = interface_residual_ML(k, omega, u, Ns, N_sub, p_size, ...
    lambda_sub, mu_sub, rho_sub, zeta_phi, zeta_psi, S0_last, S1_last)
% Analytical outgoing-wave interface check for the bottom of the last layer.

n = length(Ns);
N_last = Ns(n);
last_offset = 2 * sum(Ns(1:(n-1)));
u_last = u(last_offset + (1:2*N_last));

phi = u(p_size + (1:N_sub));
psi = u(p_size + N_sub + (1:N_sub));

A = phi(1);       % substrate potential value at the interface
B = psi(1);       % substrate shear potential value at the interface

cL = sqrt((lambda_sub + 2.0*mu_sub) / rho_sub);
cT = sqrt(mu_sub / rho_sub);

qL = outgoing_vertical_wavenumber((omega/cL)^2 - k^2, zeta_phi);
qT = outgoing_vertical_wavenumber((omega/cT)^2 - k^2, zeta_psi);

if ~(isfinite(real(qL)) && isfinite(imag(qL)) && isfinite(real(qT)) && isfinite(imag(qT)))
    bc_res = inf;
    return;
end

% Analytical derivatives at the bottom interface z = z0.
dphi  = -1i*qL*A;
ddphi = -(qL^2)*A;
dpsi  = -1i*qT*B;
ddpsi = -(qT^2)*B;

ux_sub = 1i*k*A - dpsi;
uz_sub = dphi + 1i*k*B;

tzz_sub = (lambda_sub + 2.0*mu_sub)*ddphi ...
    - lambda_sub*(k^2)*A ...
    + 2i*mu_sub*k*dpsi;

txz_sub = mu_sub*(2i*k*dphi - ddpsi - (k^2)*B);

ux_plate = u_last(N_last);
uz_plate = u_last(2*N_last);

t_plate = (S0_last + k*S1_last) * u_last;
tzz_plate = t_plate(N_last);
txz_plate = t_plate(2*N_last);

disp_jump = [ux_plate - ux_sub; uz_plate - uz_sub];
trac_jump = [tzz_plate - tzz_sub; txz_plate - txz_sub];

disp_scale = max(norm([ux_plate; uz_plate; ux_sub; uz_sub]), eps);
trac_scale = max(norm([tzz_plate; txz_plate; tzz_sub; txz_sub]), eps);

disp_res = norm(disp_jump) / disp_scale;
trac_res = norm(trac_jump) / trac_scale;

bc_res = max(disp_res, trac_res);
end

function q = outgoing_vertical_wavenumber(q2, zeta)
% Select square-root branch consistent with numerical decay along the
% rational/complex half-space map z = z0 - zeta*r, r >= 0.
q = sqrt(q2);

if imag(q*zeta) < 0
    q = -q;
elseif abs(imag(q*zeta)) <= 100*eps*max(1, abs(q*zeta))
    if real(q) < 0
        q = -q;
    end
end
end

function [Ux_layers, Uz_layers] = get_Lamb_ML_displacement_profile(U, Ns)
Ux_layers = cell(1, length(Ns));
Uz_layers = cell(1, length(Ns));

for n = 1:length(Ns)
    ind0 = 2*sum(Ns(1:(n-1)));
    Ux_layers{n} = U(ind0 + (1:Ns(n)));
    Uz_layers{n} = U(ind0 + Ns(n) + (1:Ns(n)));
end
end
