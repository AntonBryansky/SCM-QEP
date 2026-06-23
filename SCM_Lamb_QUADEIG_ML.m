%% Startup cleanup
clearvars; close all;

%% Material properties
% Desnity of material
rho = 1550.0;  % [kg/m^3]
% Elastic properties of material (stiffness matrix)
C_CFRP = [154.0, 3.7,    3.7,    0,      0,      0; % [GPa]
          3.7,   9.5,    5.2,    0,      0,      0;
          3.7,   5.2,    9.5,    0,      0,      0;
          0,     0,      0,      2.15,   0,      0;
          0,     0,      0,      0,      4.2,    0;
          0,     0,      0,      0,      0,      4.2] * 1e9;
% Thickness of one layer
h = 0.5e-3; % [m]

%% Material properties and general parameters
% Mechanical properties of layers
% Each cell/vector value corresponds to one layer from top to bottom.
densities = [rho rho rho rho]; % density, kg/m^3
Cs = {C_CFRP C_CFRP C_CFRP C_CFRP}; % stiffness matrices, Pa
orientations = [0 0 0 0]; % orientation angles of the layer, degree
thicknesses = [h h h h]; % thicknesses of the layers, m
% Differential matrices parameters, number of collocation points
Ns = [11 11 11 11];

num_of_layers = length(thicknesses); % number of layers
propagationAngle = 0; % Wave propagation direction relative to layer orientation angles, degree

freq_limit = 1e6;               % Frequency limit for plots, Hz
F_amount = 500;                 % Amount of frequency scanning points
frequency = linspace(0, freq_limit, F_amount);
frequency(1) = 1000;            % Avoid zero frequency
omegas = 2.0 * pi * frequency;  % Angular frequencies

% Wavenumber and phase velocity limits for plots and filtering
WN_limit = 4000.0;
PV_limit = 12000.0;
att_limit = 1.0;                % Maximum accepted abs(Im(k)), 1/m

%% Check QEP solver
% QEP-SCM CHANGED: the multilayer solution requires quadeig exactly as in
% the single-layer QEP implementation.
if exist('quadeig', 'file') ~= 2
    error(['quadeig was not found on the MATLAB path. ', ...
        'Install/add quadeig first, then rerun this file.']);
end

%% Transform stiffness matrix relatively to the propagation angle
beta = orientations - propagationAngle; % angle between the main direction and propagation angle, degree
c = cell(1, num_of_layers); % Set transformed matrix
for i = 1:num_of_layers
    c{1, i} = transform_stiffness_matrix(Cs{i}, beta(i));
end

%% Calculation
% Chebyshev differentiation matrices
DM1 = cell(1, num_of_layers);
DM2 = cell(1, num_of_layers);
z = cell(1, num_of_layers);
z_mm = cell(1, num_of_layers);
I = cell(1, num_of_layers);
total_h = sum(thicknesses);

for n = 1:num_of_layers
    [s, DM] = chebdif(Ns(n), 2);
    z_min = total_h/2 - sum(thicknesses(1:n));
    z_max = total_h/2 - sum(thicknesses(1:(n-1)));
    z_c = 0.5*(z_min + z_max);
    scale = 2.0/(z_max - z_min);
    DM1{n} = DM(:,:,1) * scale;
    DM2{n} = DM(:,:,2) * scale^2;
    z{n} = 0.5*(z_max - z_min)*s + z_c; % physical thickness coordinate, m
    z_mm{n} = z{n} * 1e3; % mm
    I{n} = eye(Ns(n));
end

% Mass matrix
M_layers = cell(1, num_of_layers);
for n = 1:num_of_layers
    M = get_Lamb_M_matrix(Ns(n), densities(n));
    M_layers{n} = set_Lamb_M_BC(M, Ns(n)); % Boundary/interface rows do not contain mass terms
end
M_mat = assemble_Lamb_ML_M_matrix(M_layers, Ns);

% Assemble coefficient matrices  for each prescribed wavenumber.
L0_layers = cell(1, num_of_layers);
L1_layers = cell(1, num_of_layers);
L2_layers = cell(1, num_of_layers);
S0_layers = cell(1, num_of_layers);
S1_layers = cell(1, num_of_layers);

for n = 1:num_of_layers
    [L0_layers{n}, L1_layers{n}, L2_layers{n}] = get_Lamb_L_coeffs(c{n}, I{n}, DM1{n}, DM2{n});
    [S0_layers{n}, S1_layers{n}] = get_Lamb_S_coeffs(c{n}, I{n}, DM1{n});
end

% Global multilayer QEP coefficients:
%     (A2*k^2 + A1*k + A0)U = 0,
% where A0 is later shifted by -omega^2*M.
[L0_bc, L1_bc, L2_bc] = assemble_Lamb_ML_QEP_coeffs(L0_layers, L1_layers, L2_layers, ...
                                                       S0_layers, S1_layers, I, Ns);

% Init output data storage
% Global
freqs_output = [];
cph_output = [];
att_output = [];
k_output = [];

% Symmetric modes
freqs_S_output = [];
cph_S_output = [];
att_S_output = [];
k_S_output = [];

% Antisymmetric modes
freqs_AS_output = [];
cph_AS_output = [];
att_AS_output = [];
k_AS_output = [];

% Calculation dispersion curves
timestart = tic;
f = waitbar(0, 'Calculations.', 'Name', 'Calculations...');

n_omegas = length(omegas);

for i = 1:n_omegas
    omega = omegas(i);

    % QEP-SCM CHANGED:
    %     (L0 + L1*k + L2*k^2 - omega^2*M)U = 0
    % or
    %     (A0 + A1*k + A2*k^2)U = 0
    A0 = L0_bc - (omega^2) * M_mat;
    A1 = L1_bc;
    A2 = L2_bc;

    [kvals, U] = solve_QEP_quadeig(A0, A1, A2);

    % Root filtering procedure
    for j = 1:numel(kvals)
        k = kvals(j);       % wavenumber

        % Keep only wavenumbers of propagation in the positive direction
        if real(k) <= 0
            continue;
        end
        % Keep only finite wavenumber filtering
        if ~(isfinite(real(k)) && isfinite(imag(k)))
            continue;
        end
        % Remove highly attenuated roots.
        if abs(imag(k)) > att_limit
            continue;
        end

        cph = omega / real(k);      % phase velocity, m/s
        att = imag(k);              % attenuation, 1/m
        % Keep only finite phase velocity filtering
        if cph == 0 || ~isfinite(cph)
            continue;
        end

        % Lamb mode separation
        % Note: strict symmetric/antisymmetric separation is meaningful for
        % symmetric laminate layups. For unsymmetric laminates, the global
        % output should be used without S/AS interpretation.
        u = U(:, j);
        % Rotate eigenvector because the QEP eigenvector phase is arbitrary
        [~, ind_max] = max(abs(u));
        if abs(u(ind_max)) > 0
            u = u * exp(-1i * angle(u(ind_max)));
        end

        [U11, U1n, U31, U3n] = get_Lamb_ML_border_displacements(u, Ns);
        U11 = real(U11); U1n = real(U1n);
        U31 = real(U31); U3n = real(U3n);
        er = 1e-7;

        if (abs(abs(U11 - U1n)) <= er) && (abs(abs(U31) - abs(U3n)) <= er) && (abs(U1n) >= er || (abs(U31 + U3n) <= er))
            % Symmetric mode:
            % Interface U1: displacements on the borders are equal;
            % Interface U3: modulus of displacements on the borbers are equal;
            % Displacements on the borders of U1 should be greater zero
            % or displacements on the borders of U3 should be opposite in sign
            freqs_S_output(end+1,1) = omega / (2*pi);
            cph_S_output(end+1,1) = real(cph);
            att_S_output(end+1,1) = real(att);
            k_S_output(end+1,1) = real(k);
        else
            % Otherwise the mode is antisymmetric
            freqs_AS_output(end+1,1) = omega / (2*pi);
            cph_AS_output(end+1,1) = real(cph);
            att_AS_output(end+1,1) = real(att);
            k_AS_output(end+1,1) = real(k);
        end

        % Optional plots of displacement profiles for all layers
        % [Ux_layers, Uz_layers] = get_Lamb_ML_displacement_profile(u, Ns);
        % figure(100); clf; hold on;
        % for n = 1:num_of_layers
        %     Ux_plot = real(Ux_layers{n});
        %     Uz_plot = real(Uz_layers{n});
        %     Ux_plot = Ux_plot / max(abs(Ux_plot));
        %     Uz_plot = Uz_plot / max(abs(Uz_plot));
        %     plot(Ux_plot, z_mm{n}, '-o', 'DisplayName', sprintf('Re(Ux), layer %d', n));
        %     plot(Uz_plot, z_mm{n}, '-^', 'DisplayName', sprintf('Re(Uz), layer %d', n));
        % end
        % xlabel('Normalized amplitude');
        % ylabel('Thickness coordinate z (mm)');
        % title(sprintf('f = %.2f kHz, c_p = %.1f m/s', ...
        %     omega/(2*pi)/1e3, omega/real(k)));
        % legend('Location','best'); grid on; drawnow;

        freqs_output(end+1,1) = omega / (2*pi);
        cph_output(end+1,1) = real(cph);
        att_output(end+1,1) = real(att);
        k_output(end+1,1) = real(k);
    end

    % ETA evaluation
    progr = i/n_omegas;
    timeElapsed = toc(timestart);

    if progr > 0
        etaSec = timeElapsed * (1 - progr) / progr;
    else
        etaSec = NaN;
    end

    etaText = format_eta(etaSec);
    waitbar(progr, f, sprintf('%.3f%%   ETA: %s', progr*100, etaText));
end

timeend = toc(timestart); %#ok<NASGU>
waitbar(1, f, 'Calculations completed!');
delete(f);

fprintf('QEP-SCM multilayer calculation took %.3f s\n', toc(timestart));

%% Plots
% Phase velocity versus frequency plot
figure('Color','w');
plot(freqs_output*1e-3, cph_output, '.k', 'MarkerSize', 4);
xlim([0, freq_limit*1e-3]);
ylim([0, PV_limit]);
xlabel('Frequency (kHz)');
ylabel('Phase velocity (m/s)');
grid on;

figure('Color','w'); hold on;
p1 = plot(freqs_S_output*1e-3, cph_S_output, '.r', 'MarkerSize', 4);
p2 = plot(freqs_AS_output*1e-3, cph_AS_output, '.b', 'MarkerSize', 4);
xlim([0, freq_limit*1e-3]);
ylim([0, PV_limit]);
xlabel('Frequency (kHz)');
ylabel('Phase velocity (m/s)');
legend([p1(1), p2(1)], {'SM', 'ASM'}, 'Location', 'best');
grid on;

% Wavenumber versus frequency plot
figure('Color','w');
plot(freqs_output*1e-3, k_output, '.k', 'MarkerSize', 4);
xlim([0, freq_limit*1e-3]);
ylim([0, WN_limit]);
xlabel('Frequency (kHz)');
ylabel('Wavenumber Re(k_x) (1/m)');
grid on;

figure('Color','w'); hold on;
p1 = plot(freqs_S_output*1e-3, k_S_output, '.r', 'MarkerSize', 4);
p2 = plot(freqs_AS_output*1e-3, k_AS_output, '.b', 'MarkerSize', 4);
xlim([0, freq_limit*1e-3]);
ylim([0, WN_limit]);
xlabel('Frequency (kHz)');
ylabel('Wavenumber Re(k_x) (1/m)');
legend([p1(1), p2(1)], {'SM', 'ASM'}, 'Location', 'best');
grid on;

% Attenuation versus frequency plot
figure('Color','w');
plot(freqs_output*1e-3, att_output, '.k', 'MarkerSize', 4);
xlim([0, freq_limit*1e-3]);
xlabel('Frequency (kHz)');
ylabel('Attenuation Im(k) (1/m)');
grid on;

figure('Color','w'); hold on;
p1 = plot(freqs_S_output*1e-3, att_S_output, '.r', 'MarkerSize', 4);
p2 = plot(freqs_AS_output*1e-3, att_AS_output, '.b', 'MarkerSize', 4);
xlim([0, freq_limit/1000]);
xlabel('Frequency (kHz)');
ylabel('Attenuation Im(k) (1/m)');
legend([p1(1), p2(1)], {'SM', 'ASM'}, 'Location', 'best');
grid on;

%% Funtions
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

function c = transform_stiffness_matrix(C, beta)
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
% Free-surface and interface traction rows for Lamb case.
% Row block 1: sigma_zz
% Row block 2: sigma_xz
%
% Szz = (c35*D + i*k*c13*I) Ux + (c33*D + i*k*c35*I) Uz
% Sxz = (c55*D + i*k*c15*I) Ux + (c35*D + i*k*c55*I) Uz

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

function M = get_Lamb_M_matrix(N, density)
M = -density * eye(N * 2);
end

function M = set_Lamb_M_BC(M, N)
% Rows used for external boundary conditions or interlayer continuity
% conditions must not contain inertia terms.
rows = [1, N, N+1, 2*N];
for rr = rows
    M(rr, :) = 0.0;
end
end

function [L0b, L1b, L2b] = set_coeffs_Lamb_BC(L0, L1, L2, S0, S1, N)
% Single-layer case. This part is kept for consistency and testing.
rows = [1, N, N+1, 2*N];
L0b = complex(L0);
L1b = complex(L1);
L2b = complex(L2);

for rr = rows
    L0b(rr, :) = S0(rr, :);
    L1b(rr, :) = S1(rr, :);
    L2b(rr, :) = 0.0;
end
end

function M = assemble_Lamb_ML_M_matrix(M_layers, Ns)
total_dof = 2*sum(Ns);
M = zeros(total_dof, total_dof);

for n = 1:length(Ns)
    ind = (2*sum(Ns(1:(n-1)))+1):2*sum(Ns(1:n));
    M(ind, ind) = M_layers{n};
end
end

function [A0, A1, A2] = assemble_Lamb_ML_QEP_coeffs(L0_layers, L1_layers, L2_layers, S0_layers, S1_layers, Is, Ns)
% QEP-SCM CHANGED:
% Assemble the multilayer coefficient matrices for
%
%     (A0 + A1*k + A2*k^2)U = 0.
%
% The original multilayer code assembled L(k) directly. Here each row is
% assembled in coefficient form. Bulk equations use L0, L1, L2. Free-surface
% and traction-continuity equations use S0, S1, and zero k^2 coefficient.
% Displacement-continuity equations use only the k-independent coefficient.

num_of_layers = length(Ns);
total_dof = 2*sum(Ns);
A0 = complex(zeros(total_dof, total_dof));
A1 = complex(zeros(total_dof, total_dof));
A2 = complex(zeros(total_dof, total_dof));

if num_of_layers == 1
    [A0, A1, A2] = set_coeffs_Lamb_BC(L0_layers{1}, L1_layers{1}, L2_layers{1}, ...
                                      S0_layers{1}, S1_layers{1}, Ns(1));
    return;
end

for n = 1:num_of_layers
    if n == 1 % First layer
        NL = Ns(n); NR = Ns(n + 1);

        buff_L0 = complex(L0_layers{n});
        buff_L1 = complex(L1_layers{n});
        buff_L2 = complex(L2_layers{n});
        buff_R0 = complex(zeros(2*NL, 2*NR));
        buff_R1 = complex(zeros(2*NL, 2*NR));
        buff_R2 = complex(zeros(2*NL, 2*NR));

        % External free-surface BC at the top of the first layer
        buff_L0(1, :) = S0_layers{n}(1, :);
        buff_L1(1, :) = S1_layers{n}(1, :);
        buff_L2(1, :) = 0.0;

        buff_L0(NL, :) = S0_layers{n}(NL+1, :);
        buff_L1(NL, :) = S1_layers{n}(NL+1, :);
        buff_L2(NL, :) = 0.0;

        % Traction continuity with the second layer
        buff_L0(NL+1, :) = S0_layers{n}(NL, :);
        buff_L1(NL+1, :) = S1_layers{n}(NL, :);
        buff_L2(NL+1, :) = 0.0;
        buff_R0(NL+1, :) = -S0_layers{n+1}(1, :);
        buff_R1(NL+1, :) = -S1_layers{n+1}(1, :);
        buff_R2(NL+1, :) = 0.0;

        buff_L0(2*NL, :) = S0_layers{n}(2*NL, :);
        buff_L1(2*NL, :) = S1_layers{n}(2*NL, :);
        buff_L2(2*NL, :) = 0.0;
        buff_R0(2*NL, :) = -S0_layers{n+1}(NR+1, :);
        buff_R1(2*NL, :) = -S1_layers{n+1}(NR+1, :);
        buff_R2(2*NL, :) = 0.0;

        rows = 1:2*sum(Ns(1:n));
        cols = 1:2*sum(Ns(1:(n+1)));
        A0(rows, cols) = horzcat(buff_L0, buff_R0);
        A1(rows, cols) = horzcat(buff_L1, buff_R1);
        A2(rows, cols) = horzcat(buff_L2, buff_R2);

    elseif n == num_of_layers % Last layer
        NL = Ns(n - 1); NR = Ns(n);

        buff_R0 = complex(L0_layers{n});
        buff_R1 = complex(L1_layers{n});
        buff_R2 = complex(L2_layers{n});
        buff_L0 = complex(zeros(2*NR, 2*NL));
        buff_L1 = complex(zeros(2*NR, 2*NL));
        buff_L2 = complex(zeros(2*NR, 2*NL));

        % Displacement continuity with the previous layer
        buff_L0(1, :) = get_Lamb_interface_BC(Is{n-1}(NL, :), NL, 1);
        buff_R0(1, :) = get_Lamb_interface_BC(-Is{n}(1, :), NR, 1);
        buff_L1(1, :) = 0.0; buff_R1(1, :) = 0.0;
        buff_L2(1, :) = 0.0; buff_R2(1, :) = 0.0;

        buff_L0(NR, :) = get_Lamb_interface_BC(Is{n-1}(NL, :), NL, 2);
        buff_R0(NR, :) = get_Lamb_interface_BC(-Is{n}(1, :), NR, 2);
        buff_L1(NR, :) = 0.0; buff_R1(NR, :) = 0.0;
        buff_L2(NR, :) = 0.0; buff_R2(NR, :) = 0.0;

        % External free-surface BC at the bottom of the last layer
        buff_R0(NR+1, :) = S0_layers{n}(NR, :);
        buff_R1(NR+1, :) = S1_layers{n}(NR, :);
        buff_R2(NR+1, :) = 0.0;

        buff_R0(2*NR, :) = S0_layers{n}(2*NR, :);
        buff_R1(2*NR, :) = S1_layers{n}(2*NR, :);
        buff_R2(2*NR, :) = 0.0;

        rows = (2*sum(Ns(1:(n - 1)))+1):2*sum(Ns(1:n));
        cols = (2*sum(Ns(1:(n - 2)))+1):(2*sum(Ns(1:n)));
        A0(rows, cols) = horzcat(buff_L0, buff_R0);
        A1(rows, cols) = horzcat(buff_L1, buff_R1);
        A2(rows, cols) = horzcat(buff_L2, buff_R2);

    else % Middle layers
        NL = Ns(n - 1);
        NM = Ns(n);
        NR = Ns(n + 1);

        buff_L0 = complex(zeros(2*NM, 2*NL));
        buff_L1 = complex(zeros(2*NM, 2*NL));
        buff_L2 = complex(zeros(2*NM, 2*NL));
        buff_M0 = complex(L0_layers{n});
        buff_M1 = complex(L1_layers{n});
        buff_M2 = complex(L2_layers{n});
        buff_R0 = complex(zeros(2*NM, 2*NR));
        buff_R1 = complex(zeros(2*NM, 2*NR));
        buff_R2 = complex(zeros(2*NM, 2*NR));

        % Displacement continuity with the previous layer
        buff_L0(1, :) = get_Lamb_interface_BC(Is{n-1}(NL, :), NL, 1);
        buff_M0(1, :) = get_Lamb_interface_BC(-Is{n}(1, :), NM, 1);
        buff_L1(1, :) = 0.0; buff_M1(1, :) = 0.0; buff_R1(1, :) = 0.0;
        buff_L2(1, :) = 0.0; buff_M2(1, :) = 0.0; buff_R2(1, :) = 0.0;

        buff_L0(NM, :) = get_Lamb_interface_BC(Is{n-1}(NL, :), NL, 2);
        buff_M0(NM, :) = get_Lamb_interface_BC(-Is{n}(1, :), NM, 2);
        buff_L1(NM, :) = 0.0; buff_M1(NM, :) = 0.0; buff_R1(NM, :) = 0.0;
        buff_L2(NM, :) = 0.0; buff_M2(NM, :) = 0.0; buff_R2(NM, :) = 0.0;

        % Traction continuity with the next layer
        buff_M0(NM+1, :) = S0_layers{n}(NM, :);
        buff_M1(NM+1, :) = S1_layers{n}(NM, :);
        buff_M2(NM+1, :) = 0.0;
        buff_R0(NM+1, :) = -S0_layers{n+1}(1, :);
        buff_R1(NM+1, :) = -S1_layers{n+1}(1, :);
        buff_R2(NM+1, :) = 0.0;

        buff_M0(2*NM, :) = S0_layers{n}(2*NM, :);
        buff_M1(2*NM, :) = S1_layers{n}(2*NM, :);
        buff_M2(2*NM, :) = 0.0;
        buff_R0(2*NM, :) = -S0_layers{n+1}(NR+1, :);
        buff_R1(2*NM, :) = -S1_layers{n+1}(NR+1, :);
        buff_R2(2*NM, :) = 0.0;

        rows = (2*sum(Ns(1:(n - 1)))+1):2*sum(Ns(1:n));
        cols = (2*sum(Ns(1:(n - 2)))+1):2*sum(Ns(1:(n+1)));
        A0(rows, cols) = horzcat(buff_L0, buff_M0, buff_R0);
        A1(rows, cols) = horzcat(buff_L1, buff_M1, buff_R1);
        A2(rows, cols) = horzcat(buff_L2, buff_M2, buff_R2);
    end
end
end

function I_bc = get_Lamb_interface_BC(I_n, N, Pos)
I_bc = zeros(1, 2*N);
I_bc(1, (N*(Pos-1) + 1):N*Pos) = I_n;
end

function [U11, U1n, U31, U3n] = get_Lamb_ML_border_displacements(U, Ns)
% Top surface of the first layer
U11 = U(1);
U31 = U(Ns(1) + 1);

% Bottom surface of the last layer
n = length(Ns);
ind0 = 2*sum(Ns(1:(n-1)));
U1n = U(ind0 + Ns(n));
U3n = U(ind0 + 2*Ns(n));
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

% Utilities
function etaText = format_eta(etaSec)

    if ~isfinite(etaSec) || etaSec < 0
        etaText = '--:--:--';
        return;
    end

    hours   = floor(etaSec / 3600);
    minutes = floor(mod(etaSec, 3600) / 60);
    seconds = floor(mod(etaSec, 60));

    etaText = sprintf('%02d:%02d:%02d', hours, minutes, seconds);
end
