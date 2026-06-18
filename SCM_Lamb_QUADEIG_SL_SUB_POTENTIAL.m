% Startup cleanup
clearvars; close all;

%% Material properties
% Properties of the CFRP T800/913
% Density of the CFRP T800/913
rho = 1550.0;  % [kg/m^3]
C_CFRP = [154.0, 3.7,    3.7,    0,      0,      0;  % [GPa]
          3.7,   9.5,    5.2,    0,      0,      0;
          3.7,   5.2,    9.5,    0,      0,      0;
          0,     0,      0,      2.15,   0,      0;
          0,     0,      0,      0,      4.2,    0;
          0,     0,      0,      0,      0,      4.2] * 1e9;
% Thickness of the CFRP layer
h = 2e-3; % [m]

% Properties of the aluminium substrate
% Density of the aluminium substrate
rho_sub = 2700.0; % [kg/m^3]
% Elastic properties of the aluminium substrate
C_sub = [107.5   54.59     54.59     0       0       0; % [GPa]
         54.59   107.5     54.59     0       0       0;
         54.59   54.59     107.5     0       0       0;
         0       0         0         26.45   0       0;
         0       0         0         0       26.45   0;
         0       0         0         0       0       26.45] * 1e9;

% Lame constant of the substrate
lambda_sub = C_sub(1,2);
mu_sub = C_sub(5,5);

ct_sub = sqrt(mu_sub/rho_sub);
cl_sub = sqrt((lambda_sub + 2*mu_sub)/rho_sub);

% Root filter parameters
k_abs_max = 1e7;            % Maximum value of wavenumber
att_min_tol = -1e-6;        % with exp(i*k*x-i*w*t), Im(k)>=0 is attenuating
att_max = 15e3;             % Maximum value of attenuation
qep_res_tol = 1e-6;         % residual tolerance value
interface_res_tol = 1e-5;   % interface residual tolerance value
kh_min = 0.116;             % Wavenumber value to cut off


%% Layup properties and general parameters
% Orientation angle of the layer, degrees
orientation = 0.0; % [degrees]
% Wave propagation direction angle relatively to the layer orientation
propagation_angle = 0.0; % [degrees]
% Differentiation matrices parameters, number of collocation points
N = 100; % CFRP T800/913 layer collocation points
N_sub = 100; % Aluminium substrate collocation points
% Frequency limit for plots
freq_limit = 2e6;   % [Hz]
% Phase-velocity plot limit
PV_limit = 12000.0; % [m/s]
% wavenumber limit for plots
WN_limit = 8000.0;  % rad/m for plotting y-axis
% Number of frequency samples
F_amount = 150;


%% Transform stiffness matrices relative to propagation direction
% Angle between the main direction and the plate orientations
beta = orientation - propagation_angle;
% Transform the stiffness matrix to the propagation direction
c = transform_stiffness_matrix(C_CFRP, beta);


%% Half-space map parameter.
% Real zeta gives evanescent rational-Chebyshev mapping.
% Imaginary zeta gives complex continuation for leaky radiation.
% zeta_phi controls the longitudinal potential; zeta_psi controls the shear potential.
multiplier = 20;
cases(1).name = 'evanescent';
cases(1).zeta_phi = multiplier*h;
cases(1).zeta_psi = multiplier*h;

cases(2).name = 'shear-leaky';
cases(2).zeta_phi = multiplier*h;
cases(2).zeta_psi = 1i*multiplier*h;

cases(3).name = 'fully-leaky';
cases(3).zeta_phi = 1i*multiplier*h;
cases(3).zeta_psi = 1i*multiplier*h;


%% CFRP T800/913 layer setup
% Chebyshev differentiation matrices for the CFRP layer
[s, DM] = chebdif(N, 2);
% Scaling coefficient for domain [-h/2 h/2]
scale = 2.0 / h;
DM1 = DM(:,:,1) * scale;
DM2 = DM(:,:,2) * scale^2;
% Physical thickness coordinate, m
z = (h/2) * s;
z_mm = z * 1e3;     % mm
I_N = eye(N);       % identity matrix
% Lamb polynomial coefficients for the CFRP layer
% Plate: [Ux_plate; Uz_plate]
M_plate = get_Lamb_M_matrix(N, rho);
[L0_plate, L1_plate, L2_plate] = get_Lamb_L_coeffs(c, I_N, DM1, DM2);
[S0_plate, S1_plate] = get_Lamb_S_coeffs(c, I_N, DM1);


%% Substrate layer setup and calculations

% Itialize output data storage
bc_all = [];

freq_all = [];
cph_all = [];
att_all = [];
k_all = [];
res_all = [];
case_all = [];

t0 = tic;
for case_i = 1:numel(cases)
    zeta_phi_sub = cases(case_i).zeta_phi;
    zeta_psi_sub = cases(case_i).zeta_psi;

    % Chebyshev differentiation matrices for the aluminium half-space
    z_interface = -h/2;
    [s_phi_sub, z_phi_sub, DM1_phi_sub, DM2_phi_sub] = chebdif_sub(N_sub, z_interface, zeta_phi_sub);
    [s_psi_sub, z_psi_sub, DM1_psi_sub, DM2_psi_sub] = chebdif_sub(N_sub, z_interface, zeta_psi_sub);
    z_sub_mm = z_phi_sub * 1e3;
    I_sub = eye(N_sub);

    M_sub = get_Lamb_M_matrix(N_sub, rho_sub);
    [L0_sub, L1_sub, L2_sub] = get_Lamb_potential_L_coeffs(lambda_sub, mu_sub, I_sub, DM2_phi_sub, DM2_psi_sub);
    [U0_sub, U1_sub, U2_sub] = get_Lamb_potential_U_coeffs(I_sub, DM1_phi_sub, DM1_psi_sub);
    [T0_sub, T1_sub, T2_sub] = get_Lamb_potential_T_coeffs(lambda_sub, mu_sub, I_sub, DM1_phi_sub, DM2_phi_sub, DM1_psi_sub, DM2_psi_sub);

    %% Global system assembly
    L0_global = blkdiag(L0_plate, L0_sub);
    L1_global = blkdiag(L1_plate, L1_sub);
    L2_global = blkdiag(L2_plate, L2_sub);
    M_global  = blkdiag(M_plate,  M_sub);

    % Row indices for boundary conditions
    p_i_top = [1, N+1];
    p_i_bot = [N, 2*N];
    s_i_top = [1, N_sub+1];
    s_i_bot = [N_sub, 2*N_sub];
    p_size = 2*N;
    s_size = 2*N_sub;
    total_size = p_size + s_size;

    % 1. Top free surface conditions of the T800/913 CFRP layer
    for j = 1:2
        L0_global(p_i_top(j), :) = 0.0;
        L1_global(p_i_top(j), :) = 0.0;
        L2_global(p_i_top(j), :) = 0.0;
        M_global(p_i_top(j), :)  = 0.0;

        L0_global(p_i_top(j), 1:p_size) = S0_plate(p_i_top(j), :);
        L1_global(p_i_top(j), 1:p_size) = S1_plate(p_i_top(j), :);
    end

    % 2. Interface traction continuity at z = -h/2:
    for j = 1:2
        L0_global(p_i_bot(j), :) = 0.0;
        L1_global(p_i_bot(j), :) = 0.0;
        L2_global(p_i_bot(j), :) = 0.0;
        M_global(p_i_bot(j), :)  = 0.0;

        L0_global(p_i_bot(j), 1:p_size) = S0_plate(p_i_bot(j), :);
        L1_global(p_i_bot(j), 1:p_size) = S1_plate(p_i_bot(j), :);

        L0_global(p_i_bot(j), p_size + (1:s_size)) = -T0_sub(s_i_top(j), :);
        L1_global(p_i_bot(j), p_size + (1:s_size)) = -T1_sub(s_i_top(j), :);
        L2_global(p_i_bot(j), p_size + (1:s_size)) = -T2_sub(s_i_top(j), :);
    end

    % 3. Interface displacement continuity at z = -h/2:
    for j = 1:2
        L0_global(p_size + s_i_top(j), :) = 0.0;
        L1_global(p_size + s_i_top(j), :) = 0.0;
        L2_global(p_size + s_i_top(j), :) = 0.0;
        M_global(p_size + s_i_top(j), :)  = 0.0;

        L0_global(p_size + s_i_top(j), p_i_bot(j)) = 1.0;

        L0_global(p_size + s_i_top(j), p_size + (1:s_size)) = -U0_sub(s_i_top(j), :);
        L1_global(p_size + s_i_top(j), p_size + (1:s_size)) = -U1_sub(s_i_top(j), :);
        L2_global(p_size + s_i_top(j), p_size + (1:s_size)) = -U2_sub(s_i_top(j), :);
    end

    % 4. Substrate condition at infinity
    for j = 1:2
        L0_global(p_size + s_i_bot(j), :) = 0.0;
        L1_global(p_size + s_i_bot(j), :) = 0.0;
        L2_global(p_size + s_i_bot(j), :) = 0.0;
        M_global(p_size + s_i_bot(j), :)  = 0.0;

        L0_global(p_size + s_i_bot(j), p_size + s_i_bot(j)) = 1.0;
    end

    % Frequency sweep
    freqs = linspace(1.0, freq_limit, F_amount);
    omegas = 2.0 * pi * freqs;

    freq_plot = [];
    cph_plot = [];
    att_plot = [];
    k_plot = [];
    res_plot = [];
    U_modes = [];

    bc_res_plot = [];

    fprintf('Starting Lamb QEP calculation for CFRP layer on aluminium substrate...\n');
    for i = 1:numel(omegas)
        omega = omegas(i);

        if mod(i-1, 10) == 0
            pct = 100.0 * (i-1) / numel(omegas);
            fprintf('  %.1f %%\n', pct);
        end

        A0 = L0_global - (omega^2) * M_global;
        A1 = L1_global;
        A2 = L2_global;

        [kvals, U] = solve_QEP_quadeig(A0, A1, A2);

        for m = 1:numel(kvals)
            k = kvals(m);
            
            % Keep only wavenumbers of propagation in the positive direction
            if real(k) <= 0
                continue;
            end

            % Keep only finite wavenumber filtering
            if ~(isfinite(real(k)) && isfinite(imag(k)))
                continue;
            end
            
            % Filtering relatively to imaginary wavenumber
            alpha = imag(k);
            % Keep only roots with non-negative attenuation
            if alpha < att_min_tol
                continue;
            end
            % Remove highly attenuated roots.
            if alpha > att_max
                continue
            end

            if abs(k) > k_abs_max
                continue;
            end

            cph = omega / real(k);

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
                kh_val = real(k) * h;
                if ~(isfinite(kh_val) && kh_val >= kh_min)
                    continue;
                end
            end

            if ~(isfinite(cph) && cph > 0)
                continue;
            end

            u = U(:, m);
            if any(~isfinite(u))
                continue;
            end

            res = qep_residual(A0, A1, A2, k, u);
            if res > qep_res_tol
                continue;
            end

            

            bc_res = interface_residual( ...
                k, omega, u, N, N_sub, p_size, ...
                lambda_sub, mu_sub, rho_sub, ...
                zeta_phi_sub, zeta_psi_sub, S0_plate, S1_plate);
            if ~(isfinite(bc_res) && bc_res <= interface_res_tol)
                continue;
            end

            % PLoting profiles
            u_layer = u(1:p_size);

            % Rotate eigenvector because the QEP eigenvector phase is arbitrary.
            [~, ind_max] = max(abs(u_layer));
            if abs(u_layer(ind_max)) > 0
                u_layer = u_layer * exp(-1i * angle(u_layer(ind_max)));
            end

            Ux = u_layer(1:N);
            Uz = u_layer(N+1:2*N);

            % Normalization for simple profile inspection.
            Ux_norm = max(abs(real(Ux)));
            Uz_norm = max(abs(real(Uz)));

            if Ux_norm > 0
                Ux = Ux / Ux_norm;
            end
            if Uz_norm > 0
                Uz = Uz / Uz_norm;
            end

            freq_plot(end+1,1) = omega / (2*pi); %#ok<AGROW>
            cph_plot(end+1,1) = real(cph); %#ok<AGROW>
            att_plot(end+1,1) = imag(k); %#ok<AGROW>
            k_plot(end+1,1) = k; %#ok<AGROW>
            res_plot(end+1,1) = res; %#ok<AGROW>
            bc_res_plot(end+1,1) = bc_res; %#ok<AGROW>
            U_modes(:,end+1) = u; %#ok<AGROW>
        end
    end
    freq_all = [freq_all; freq_plot];
    cph_all = [cph_all; cph_plot];
    att_all = [att_all; att_plot];
    k_all = [k_all; k_plot];
    res_all = [res_all; res_plot];
    bc_all   = [bc_all;   bc_res_plot];
    case_all = [case_all; case_i * ones(size(freq_plot))];

    fprintf('Case %s accepted roots: %d\n', cases(case_i).name, numel(k_plot));
end

fprintf('  100.0 %%\n');
fprintf('QEP calculation took %.3f s\n', toc(t0));
fprintf('Accepted roots: %d\n', numel(k_all));

%% Plots
k_real   = real(k_all);           % [1/m]
alpha    = imag(k_all);           % [1/m]

% Bulk-wave threshold lines in k-f space:
% cph = omega/k  =>  k = omega/c
freq_line = linspace(0, freq_limit, 500).';
omega_line = 2*pi*freq_line*1e3;

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

% Phase velocity versus frequency plot
figure('Color','w');
plot(freq_all*1e-3, cph_all, '.r', 'MarkerSize', 6);
yline(ct_sub, '--k', 'c_T Al');
yline(cl_sub, '--k', 'c_L Al');
xlim([0, freq_limit*1e-3]);
ylim([0, PV_limit]);
xlabel('Frequency (kHz)');
ylabel('Phase velocity (m/s)');
grid on;


% Attenuation versus frequency plot
figure('Color','w'); 
hold on;

idx = case_all == 1;
plot(freq_all(idx)*1e-3, att_all(idx), '.r', ...
    'MarkerSize', 6, ...
    'DisplayName', 'evanescent');

idx = case_all == 2;
plot(freq_all(idx)*1e-3, att_all(idx), '.b', ...
    'MarkerSize', 6, ...
    'DisplayName', 'shear-leaky');

idx = case_all == 3;
plot(freq_all(idx)*1e-3, att_all(idx), '.k', ...
    'MarkerSize', 6, ...
    'DisplayName', 'fully-leaky');

xlim([0, freq_limit*1e-3]);
ylim([0, att_max]);

xlabel('Frequency (kHz)');
ylabel('Attenuation Im(k_x) (1/m)');

legend('Location','best');
grid on;
box on;


%% Functions
function [kvals, U] = solve_QEP_quadeig(A0, A1, A2)
% This function solves the quadratic eigenvalue problem (QEP):
%
%     (A2*k^2 + A1*k + A0) U = 0
%
% where:
%     k  - unknown wavenumber
%     U  - corresponding eigenvector / displacement-type solution
%
% The problem is solved using the external function quadeig.
%
% To improve numerical conditioning, the original wavenumber k is scaled as
%
%     k = k0 * lambda
%
% Therefore, the QEP is rewritten in terms of lambda:
%
%     (B2*lambda^2 + B1*lambda + B0) U = 0
%
% where:
%
%     B0 = A0
%     B1 = A1*k0
%     B2 = A2*k0^2
%
% After solving for lambda, the physical wavenumber is recovered as:
%
%     k = k0 * lambda

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
%
% Lamb/sagittal block:
%     L11 = c55*D2 - c11*I*k^2 + 2*i*c15*k*D
%     L13 = c35*D2 - c15*I*k^2 + i*(c13+c55)*k*D
%     L31 = L13
%     L33 = c33*D2 - c55*I*k^2 + 2*i*c35*k*D

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
%
% sigma_zz = (c35*D + i*k*c13*I) Ux + (c33*D + i*k*c35*I) Uz
% sigma_xz = (c55*D + i*k*c15*I) Ux + (c35*D + i*k*c55*I) Uz

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
%
% phi equation:
%     (lambda + 2*mu)*(D2 - k^2*I)*phi + rho*omega^2*phi = 0
% psi equation:
%     mu*(D2 - k^2*I)*psi + rho*omega^2*psi = 0

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
%
% sigma_zz = ((lambda + 2*mu)*D2 - lambda*k^2*I)*phi
%            + 2*i*mu*k*D(psi)
%
% sigma_xz = mu*(2*i*k*D(phi) - D2(psi) - k^2*psi)

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

function r = qep_residual(A0, A1, A2, k, u)
num = norm((A2*k^2 + A1*k + A0) * u);
den = (norm(A2)*abs(k)^2 + norm(A1)*abs(k) + norm(A0)) * norm(u);
if den == 0
    r = inf;
else
    r = num / den;
end
end

function bc_res = interface_residual(k, omega, u, N, N_sub, p_size, ...
    lambda_sub, mu_sub, rho_sub, zeta_phi, zeta_psi, S0_plate, S1_plate)
% Analytical outgoing-wave interface check.
%
% The QEP enforces the interface with the Chebyshev derivative matrices.
% This additional test reconstructs the substrate field at the interface
% from the analytical outgoing half-space potentials
%
%   phi(z) = A exp[-i*qL*(z-z0)],   psi(z) = B exp[-i*qT*(z-z0)],
%
% where z0 = -h/2 is the interface and z <= z0 is the substrate.  The square
% root branches qL and qT are selected so that the field decays along the
% rational/complex map used in the SCM substrate discretisation.  Roots that
% satisfy the matrix equation but do not satisfy this analytical outgoing
% interface form are rejected as bulk/non-converged/spurious roots.

u_plate = u(1:p_size);
phi = u(p_size + (1:N_sub));
psi = u(p_size + N_sub + (1:N_sub));

A = phi(1);       % substrate potential value at the interface
B = psi(1);       % substrate shear potential value at the interface

% Bulk wavenumbers of the isotropic substrate.
cL = sqrt((lambda_sub + 2.0*mu_sub) / rho_sub);
cT = sqrt(mu_sub / rho_sub);

qL = outgoing_vertical_wavenumber((omega/cL)^2 - k^2, zeta_phi);
qT = outgoing_vertical_wavenumber((omega/cT)^2 - k^2, zeta_psi);

if ~(isfinite(real(qL)) && isfinite(imag(qL)) && isfinite(real(qT)) && isfinite(imag(qT)))
    bc_res = inf;
    return;
end

% Analytical derivatives at the bottom interface z = z0.
% The lower half-space outgoing form is exp[-i*q*(z-z0)].
dphi  = -1i*qL*A;
ddphi = -(qL^2)*A;
dpsi  = -1i*qT*B;
ddpsi = -(qT^2)*B;

% Substrate displacement reconstructed analytically from potentials.
ux_sub = 1i*k*A - dpsi;
uz_sub = dphi + 1i*k*B;

% Substrate tractions reconstructed analytically from potentials.
tzz_sub = (lambda_sub + 2.0*mu_sub)*ddphi ...
    - lambda_sub*(k^2)*A ...
    + 2i*mu_sub*k*dpsi;

txz_sub = mu_sub*(2i*k*dphi - ddpsi - (k^2)*B);

% Plate displacement and traction at the bottom interface.
ux_plate = u_plate(N);
uz_plate = u_plate(2*N);

t_plate = (S0_plate + k*S1_plate) * u_plate;
tzz_plate = t_plate(N);
txz_plate = t_plate(2*N);

disp_jump = [ux_plate - ux_sub; uz_plate - uz_sub];
trac_jump = [tzz_plate - tzz_sub; txz_plate - txz_sub];

disp_scale = max(norm([ux_plate; uz_plate; ux_sub; uz_sub]), eps);
trac_scale = max(norm([tzz_plate; txz_plate; tzz_sub; txz_sub]), eps);

disp_res = norm(disp_jump) / disp_scale;
trac_res = norm(trac_jump) / trac_scale;

bc_res = max(disp_res, trac_res);
end

function q = outgoing_vertical_wavenumber(q2, zeta)
% Select the square-root branch consistent with numerical decay along the
% rational/complex half-space map z = z0 - zeta*r, r >= 0.
%
% For the lower half-space outgoing potential exp[-i*q*(z-z0)], the field
% sampled along the mapped path behaves as exp[i*q*zeta*r].  Numerical
% decay requires imag(q*zeta) > 0.  When the test is almost neutral, use the
% branch with positive real(q) as a stable tie-breaker.

q = sqrt(q2);

if imag(q*zeta) < 0
    q = -q;
elseif abs(imag(q*zeta)) <= 100*eps*max(1, abs(q*zeta))
    if real(q) < 0
        q = -q;
    end
end
end