% Startup cleanup
clearvars; close all;

%% Material properties
% Desnity of the CFRP T800/913
rho = 1550.0;  % [kg/m^3]
% Elastic properties of the CFRP T800/913 (stiffness matrix)
C_CFRP = [154.0, 3.7,    3.7,    0,      0,      0; % [GPa]
          3.7,   9.5,    5.2,    0,      0,      0;
          3.7,   5.2,    9.5,    0,      0,      0;
          0,     0,      0,      2.15,   0,      0;
          0,     0,      0,      0,      4.2,    0;
          0,     0,      0,      0,      0,      4.2] * 1e9;
% Thickness of plate 
h = 2e-3; % [m]


%% Plate properties and general parameters
% Mechanical properties of the layer
density = rho;
% Orientation angle of the layer, degrees
orientation = 0.0;
% Wave propagation directionrealtive to layer orietation
propagation_angle = 0.0;

% Differentiation matrices parameter, number of collocation points
N = 21;
% Frequency limit for plot
freq_limit = 1e6;               % Hz
% Wavenumber limit 
WN_limit = 4000.0;              % rad/m for plotting y-axis
PV_limit = 12000.0;             % Limit of the phase velocity for plots, used for plot
% Amount of points per frequency axis
F_amount = 500;                 % Amount of frequency scanning point
frequency = linspace(0, freq_limit, F_amount);
frequency(1) = 1000;            % Avoid zero values
omegas = 2.0 * pi * frequency;  % Angular frequencies


%% Transform stiffness matrix relatively to the propagation direction (angle)
% Angle between the main direction and the plate orientations
beta = orientation - propagation_angle;
% Transform the stiffness matrix to the propagation direction
c = transform_stiffness_matrix(C_CFRP, beta);


%% Calculations
% Chebyshev differentiation matrices
% Calculation of the differentiation matrices for domain [-1 1]
[s, DM] = chebdif(N, 2);
% Scaling coefficient for domain [-h/2 h/2]
scale = 2.0 / h;
DM1 = DM(:,:,1) * scale;
DM2 = DM(:,:,2) * scale^2;

% Physical thickness coordinate, m
z = (h/2) * s; 
z_mm = z * 1e3; % mm
I_N = eye(N);

% Mass matrix and polynomial coefficients

n_dof = 3 * N;
M_mat = get_M_matrix(N, density);           % M matrix does not depend on the wavenumber
M_mat = set_M_BC(M_mat, N);                 % Setting boundary conditions for M matrix

% L matrices and setting boundary conditions
[L0, L1, L2] = get_L_coeffs(c, I_N, DM1, DM2);
[S0, S1] = get_S_coeffs(c, I_N, DM1);
[L0_bc, L1_bc, L2_bc] = set_coeffs_BC(L0, L1, L2, S0, S1, N);

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

% Calculation of dispersion curves
timestart = tic;                            % initialize timer
f = waitbar(0, 'Calculations.', 'Name', 'Calculations...');

n_omegas = length(omegas);

for i = 1:n_omegas
    omega = omegas(i);
    
    % (A0 + A1*k + A2*k^2)U = 0
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
        if abs(imag(k)) > 1
            continue;
        end

        cph = omega / real(k);      % phase velocity, m/s
        att = imag(k);              % attenuation, 1/m
        % Keep only finite phase velocity filtering
        if cph == 0 || ~isfinite(cph)
            continue;
        end
        
        % Mode separation
        % Displacement profile components
        u = U(:, j);
        if numel(u) ~= n_dof
            continue;
        end

        % Rotate eigenvector because the QEP eigenvector phase is arbitrary
        [~, ind_max] = max(abs(u));
        if abs(u(ind_max)) > 0
            u = u * exp(-1i * angle(u(ind_max)));
        end

        
        Ux = u(1:N);
        Uy = u(N+1:2*N);
        Uz = u(2*N+1:3*N);

        U11 = real(Ux(1)); U1n = real(Ux(end));
        U21 = real(Uy(1)); U2n = real(Uy(end));
        U31 = real(Uz(1)); U3n = real(Uz(end));
        er = 1e-7;

        % General parity separation from the coupled SCM_SL.m case.
        % Symmetric mode: Ux and Uy have the same sign at both faces,
        % while Uz is opposite-sign or close to zero at one face.
        if (U11*U1n >= -er) && (U21*U2n >= -er) && (U31*U3n <= er)
            freqs_S_output(end+1,1) = omega / (2*pi);
            cph_S_output(end+1,1) = real(cph);
            att_S_output(end+1,1) = real(att);
            k_S_output(end+1,1) = real(k);
        else
            freqs_AS_output(end+1,1) = omega / (2*pi);
            cph_AS_output(end+1,1) = real(cph);
            att_AS_output(end+1,1) = real(att);
            k_AS_output(end+1,1) = real(k);
        end

        % Optional plots of displacement profiles
        % Ux = Ux / max(abs(real(Ux)));
        % Uy = Uy / max(abs(real(Uy)));
        % Uz = Uz / max(abs(real(Uz)));
        % figure(100);
        % clf;
        % plot(real(Ux), z_mm, '-o', 'DisplayName','Re(Ux)');
        % hold on;
        % plot(real(Uy), z_mm, '-s', 'DisplayName','Re(Uy)');
        % plot(real(Uz), z_mm, '-^', 'DisplayName','Re(Uz)');
        % xlabel('Normalized amplitude');
        % ylabel('Thickness coordinate z (mm)');
        % title(sprintf('f = %.2f kHz, c_p = %.1f m/s', ...
        %     omega/(2*pi)/1e3, omega/real(k)));
        % legend('Location','best');
        % grid on;
                % drawnow;


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

timeend = toc(timestart);
waitbar(1, f, 'Calculations completed!');
delete(f);

fprintf('QEP calculation took %.3f s\n', timeend);


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


% Wavenumenber versus frequency plot
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
% Free-surface traction rows for Lamb case.
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
% Boundary rows for [Ux;Uz] are:
%   row 1      : top sigma_zz = 0
%   row N      : bottom sigma_zz = 0
%   row N+1    : top sigma_xz = 0
%   row 2N     : bottom sigma_xz = 0
rows = [1, N, N+1, 2*N];
for rr = rows
    M(rr, :) = 0.0;
end
end

function [L0b, L1b, L2b] = set_coeffs_Lamb_BC(L0, L1, L2, S0, S1, N)
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

function [L0, L1, L2] = get_L_coeffs(c, I, D, D2)
% Full coupled coefficient matrices from SCM_SL.m.
% Displacement order: [Ux; Uy; Uz]
%
% The original matrix L(k) is written as:
%
%     L(k) = L0 + L1*k + L2*k^2
%
% Therefore, for fixed omega the QEP is:
%
%     (L0 - omega^2*M + L1*k + L2*k^2) U = 0

L11_0 = c(5,5) * D2;
L12_0 = c(4,5) * D2;
L13_0 = c(3,5) * D2;
L21_0 = L12_0;
L22_0 = c(4,4) * D2;
L23_0 = c(3,4) * D2;
L31_0 = L13_0;
L32_0 = L23_0;
L33_0 = c(3,3) * D2;

L0 = [L11_0, L12_0, L13_0;
      L21_0, L22_0, L23_0;
      L31_0, L32_0, L33_0];

L11_1 = 2i * c(1,5) * D;
L12_1 = 1i * (c(1,4) + c(5,6)) * D;
L13_1 = 1i * (c(1,3) + c(5,5)) * D;
L21_1 = L12_1;
L22_1 = 2i * c(4,6) * D;
L23_1 = 1i * (c(3,6) + c(4,5)) * D;
L31_1 = L13_1;
L32_1 = L23_1;
L33_1 = 2i * c(3,5) * D;

L1 = [L11_1, L12_1, L13_1;
      L21_1, L22_1, L23_1;
      L31_1, L32_1, L33_1];

L11_2 = -c(1,1) * I;
L12_2 = -c(1,6) * I;
L13_2 = -c(1,5) * I;
L21_2 = L12_2;
L22_2 = -c(6,6) * I;
L23_2 = -c(5,6) * I;
L31_2 = L13_2;
L32_2 = L23_2;
L33_2 = -c(5,5) * I;

L2 = [L11_2, L12_2, L13_2;
      L21_2, L22_2, L23_2;
      L31_2, L32_2, L33_2];
end

function [S0, S1] = get_S_coeffs(c, I, D)
% Full coupled free-surface traction coefficient matrices from SCM_SL.m.
% The original traction matrix S(k) is written as:
%
%     S(k) = S0 + S1*k

S1_0 = c(3,5) * D;
S2_0 = c(3,4) * D;
S3_0 = c(3,3) * D;
S4_0 = c(4,5) * D;
S5_0 = c(4,4) * D;
S6_0 = c(3,4) * D;
S7_0 = c(5,5) * D;
S8_0 = c(4,5) * D;
S9_0 = c(3,5) * D;

S0 = [S1_0, S2_0, S3_0;
      S4_0, S5_0, S6_0;
      S7_0, S8_0, S9_0];

S1_1 = 1i * c(1,3) * I;
S2_1 = 1i * c(3,6) * I;
S3_1 = 1i * c(3,5) * I;
S4_1 = 1i * c(1,4) * I;
S5_1 = 1i * c(4,6) * I;
S6_1 = 1i * c(4,5) * I;
S7_1 = 1i * c(1,5) * I;
S8_1 = 1i * c(5,6) * I;
S9_1 = 1i * c(5,5) * I;

S1 = [S1_1, S2_1, S3_1;
      S4_1, S5_1, S6_1;
      S7_1, S8_1, S9_1];
end

function M = get_M_matrix(N, density)
M = -density * eye(N * 3);
end

function M = set_M_BC(M, N)
rows = [1, N, N+1, 2*N, 2*N+1, 3*N];
for rr = rows
    M(rr, :) = 0.0;
end
end

function [L0b, L1b, L2b] = set_coeffs_BC(L0, L1, L2, S0, S1, N)
rows = [1, N, N+1, 2*N, 2*N+1, 3*N];
L0b = complex(L0);
L1b = complex(L1);
L2b = complex(L2);

for rr = rows
    L0b(rr, :) = S0(rr, :);
    L1b(rr, :) = S1(rr, :);
    L2b(rr, :) = 0.0;
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
