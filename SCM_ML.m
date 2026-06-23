%% Startup cleanup
clearvars; close all;
% Coupled multilayer SCM calculation.
% The displacement vector keeps all three components: U = [U1; U2; U3].

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
% Thickness of plate 
h = 0.5e-3; % [m]

%% Material properties and general parameters
% Mechanical properties of layers
densities = [rho rho rho rho]; % density, kg/m^3
Cs = {C_CFRP C_CFRP C_CFRP C_CFRP}; % stiffness matrices, Pa
orientations = [0 0 0 0]; % orientation angles of the layer, degree
thicknesses = [h h h h]; % thicknesses of the layers, m
% Differential matrices parameters, number of collocation points
Ns = [11 11 11 11];

num_of_layers = length(thicknesses); % number of layers
propagationAngle = 0; % Wave propagation direction relative to layer orientation angles, degree

freq_limit = 1e6; % Frequency limit for plots, Hz
% Set range of wavenumbers
WN_limit = 4000;
WN_amount = 200;
wavenumber = 0:(WN_limit/WN_amount):WN_limit; wavenumber(1) = 1e-6;

PV_limit = 12000.0;

%% Transform stiffness matrix relatively to the propagation angle
beta = orientations - propagationAngle; % angle between the main direction and propagation angle, usually equals to zero, degree
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

% Independent matrices
M_layers = cell(1, num_of_layers);
for n = 1:num_of_layers
    M = get_M_matrix(Ns(n), densities(n)); % It is possible to exclude M matrix calculation from cycle.
    M_layers{n} = set_M_BC(M, Ns(n)); % Boundary conditions for M matrix
end
M = assemble_ML_M_matrix(M_layers, Ns);

% Init output phase velocities
total_N = sum(Ns);
total_dof = 3*total_N;
freqs = NaN(length(wavenumber), total_dof);
PV = NaN(length(wavenumber), total_dof);
PV_S = NaN(length(wavenumber), total_dof);
PV_AS = NaN(length(wavenumber), total_dof);
k_plot =  NaN(length(wavenumber), total_dof);

% Calculation dispersion curves
timestart = tic;

n_wavenumber = length(wavenumber);
f = waitbar(0, 'Calculations.', 'Name', 'Calculations...');

for i = 1:n_wavenumber
    % calculating L and S matrices
    k = wavenumber(i);
    L_layers = cell(1, num_of_layers);
    S_layers = cell(1, num_of_layers);

    for n = 1:num_of_layers
        L_layers{n} = get_L_matrix(c{n}, k, I{n}, DM1{n}, DM2{n});
        S_layers{n} = get_S_matrix(c{n}, k, I{n}, DM1{n});
    end

    % Setting boundary and interface conditions
    % Boundary conditions for L matrix
    L = assemble_ML_L_matrix(L_layers, S_layers, I, Ns);

    % Balancing
    [T1, T2] = balance2(L, M);
    L_bal = T1*L*T2;
    M_bal = T1*M*T2;
    [U, w2] = eig(L_bal, M_bal, 'qz');
    U = T2 * U;
    w2 = diag(w2);              % squared angular frequency
    [w2, inds] = sort(w2);
    U = U(:, inds);             
    w = real(sqrt(w2));         % angular frequency, 1/s
    Vp = w./k;                  % phase velocity, m/s
    Vp = Vp';
    fs = w./(2*pi);
    freqs(i, :) = fs/1000;
    PV(i, :) = Vp;
    k_plot(i, :) = k;

    % Mode separation
    for j = 1:total_dof
        % The eigenvectors are generally complex because L contains 1i*k terms.
        % Therefore, do not compare signs directly. Instead, compare parity
        % residuals of the top and bottom surface displacements.
        if Vp(j) ~= 0 && isfinite(Vp(j)) % Skip values equal to zero or Inf
            Uj = U(:, j);

            % Normalize the eigenvector. The absolute amplitude of an
            % eigenvector is arbitrary, but the relative parity is meaningful.
            norm_Uj = max(abs(Uj));
            if norm_Uj > 0
                Uj = Uj ./ norm_Uj;
            end

            [U11, U1n, U21, U2n, U31, U3n] = get_ML_border_displacements(Uj, Ns);

            % Symmetric Lamb-type parity for the coupled 3-component case:
            % U1 and U2 are even, U3 is odd with respect to the mid-plane.
            sym_coeff = norm([U11 - U1n, U21 - U2n, U31 + U3n]);

            % Antisymmetric Lamb-type parity for the coupled 3-component case:
            % U1 and U2 are odd, U3 is even with respect to the mid-plane.
            asym_coeff = norm([U11 + U1n, U21 + U2n, U31 - U3n]);

            if sym_coeff <= asym_coeff
                PV_S(i, j) = Vp(j);
            else
                PV_AS(i, j) = Vp(j);
            end
        end
    end


    PV(i, :) = Vp;

    % ETA evaluation
    progr = i/n_wavenumber;
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

fprintf('SCM calculation took %.3f s\n', timeend);

%% Plot

figure('Color','w');
p1 = plot(freqs, PV_S, '.r');
p1 = p1(1);
ylim([0 1.2e4]);
xlim([0 1000]);
hold on;
p2 = plot(freqs, PV_AS, '.b');
p2 = p2(1);
ylim([0 PV_limit]);
xlim([0 freq_limit/1000]);
xlabel('Frequency, kHz');
ylabel('Phase velocity, m/s');
legend([p1, p2], {'SM', 'ASM'});



figure('Color','w');
plot(freqs, PV, '.r');
ylim([0 PV_limit]);
xlim([0 freq_limit/1000]);
xlabel('Frequency, kHz');
ylabel('Phase velocity, m/s');
title('All coupled modes');



% Wavenumber-frequency plot for conventional fixed-k solution

Kmat = repmat(wavenumber(:), 1, total_dof);   % each row has the same input k


% Use your already separated symmetric/antisymmetric masks
K_S  = Kmat;
K_AS = Kmat;

K_S(~isfinite(PV_S))   = NaN;
K_AS(~isfinite(PV_AS)) = NaN;

figure('Color','w'); hold on;

p1 = plot(freqs, K_S,  '.r', 'MarkerSize', 4);
p2 = plot(freqs, K_AS, '.b', 'MarkerSize', 4);

xlabel('Frequency (kHz)');
ylabel('Wavenumber Re(k_x) (1/m)');

xlim([0, freq_limit/1000]);
ylim([0, WN_limit]);

legend([p1(1), p2(1)], {'SM', 'ASM'}, 'Location', 'best');
grid on;



K_all = Kmat;
K_all(~isfinite(PV)) = NaN;

figure('Color','w'); hold on;
plot(freqs, K_all, '.r', 'MarkerSize', 4);

xlabel('Frequency (kHz)');
ylabel('Wavenumber Re(k_x) (1/m)');

xlim([0, freq_limit/1000]);
ylim([0, WN_limit]);
grid on;


%% Funtions
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
c = zeros(6,6);
s = sind(beta);
g = cosd(beta);
c(1,1) = C(1,1)*g^4+C(2,2)*s^4+2*(C(1,2)+2*C(6,6))*s^2*g^2; %#ok<*SAGROW>
c(1,2) = (C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2+C(1,2);
c(1,3) = C(1,3)*g^2+C(2,3)*s^2;
c(1,6) = (C(1,2)+2*C(6,6)-C(1,1))*s*g^3+(C(2,2)-C(1,2)-2*C(6,6))*g*s^3;
c(2,2) = C(1,1)*s^4+C(2,2)*g^4+2*(C(1,2)+2*C(6,6))*s^2*g^2;
c(2,3) = C(2,3)*g^2+C(1,3)*s^2;
c(2,6) = (C(1,2)+2*C(6,6)-C(1,1))*g*s^3+(C(2,2)-C(1,2)-2*C(6,6))*s*g^3;
c(3,3) = C(3,3);
c(3,6) = (C(2,3)-C(1,3))*s*g;
c(4,4) = C(4,4)*g^2+C(5,5)*s^2;
c(4,5) = (C(4,4)-C(5,5))*s*g;
c(5,5) = C(5,5)*g^2+C(4,4)*s^2;
c(6,6) = C(6,6)+(C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2;
end

function L = get_L_matrix(c, k, I, D, D2)
L11 = c(5,5)*D2 - c(1,1)*I*k^2 + 2*c(1,5)*k*1i*D;
L12 = c(4,5)*D2 - c(1,6)*I*k^2 + (c(1,4) + c(5,6))*k*1i*D;
L13 = c(3,5)*D2 - c(1,5)*I*k^2 + (c(1,3) + c(5,5))*k*1i*D;
L21 = L12;
L22 = c(4,4)*D2 - c(6,6)*I*k^2 + 2*c(4,6)*k*1i*D;
L23 = c(3,4)*D2 - c(5,6)*I*k^2 + (c(3,6) + c(4,5))*k*1i*D;
L31 = L13;
L32 = L23;
L33 = c(3,3)*D2 - c(5,5)*I*k^2 + 2*c(3,5)*k*1i*D;
L = [L11 L12 L13;...
     L21 L22 L23;...
     L31 L32 L33];
end

function S = get_S_matrix(c, k, I, D)
S1 = c(1,3)*k*1i*I + c(3,5)*D;
S2 = c(3,6)*k*1i*I + c(3,4)*D;
S3 = c(3,5)*k*1i*I + c(3,3)*D;
S4 = c(1,4)*k*1i*I + c(4,5)*D;
S5 = c(4,6)*k*1i*I + c(4,4)*D;
S6 = c(4,5)*k*1i*I + c(3,4)*D;
S7 = c(1,5)*k*1i*I + c(5,5)*D;
S8 = c(5,6)*k*1i*I + c(4,5)*D;
S9 = c(5,5)*k*1i*I + c(3,5)*D;
S = [S1 S2 S3;...
     S4 S5 S6;...
     S7 S8 S9];
end

function M = get_M_matrix(N, p)
M = eye(N*3).*(-p);
end

function M = set_M_BC(M, N)
% Boundary rows: top and bottom collocation points for U1, U2 and U3.
rows = [1, N, N + 1, 2*N, 2*N + 1, 3*N];
for i = 1:length(rows)
    M(rows(i), rows(i)) = 0;
end
end

function M = assemble_ML_M_matrix(M_layers, Ns)
total_dof = 3*sum(Ns);
M = zeros(total_dof, total_dof);

for n = 1:length(Ns)
    ind = get_layer_ind(Ns, n);
    M(ind, ind) = M_layers{n};
end
end

function L = assemble_ML_L_matrix(L_layers, S_layers, Is, Ns)
num_of_layers = length(Ns);
total_dof = 3*sum(Ns);
L = zeros(total_dof, total_dof);

% First, insert all differential equations into the block diagonal matrix.
for n = 1:num_of_layers
    ind = get_layer_ind(Ns, n);
    L(ind, ind) = L_layers{n};
end

if num_of_layers == 1
    N = Ns(1);
    ind = get_layer_ind(Ns, 1);
    rows_top = [1, N + 1, 2*N + 1];
    rows_bot = [N, 2*N, 3*N];

    L(ind(rows_top), :) = 0;
    L(ind(rows_top), ind) = S_layers{1}(rows_top, :);

    L(ind(rows_bot), :) = 0;
    L(ind(rows_bot), ind) = S_layers{1}(rows_bot, :);
    return;
end

for n = 1:num_of_layers
    if n == 1 % First layer
        NL = Ns(n);
        NR = Ns(n + 1);
        ind_L = get_layer_ind(Ns, n);
        ind_R = get_layer_ind(Ns, n + 1);
        rows_top_L = [1, NL + 1, 2*NL + 1];
        rows_bot_L = [NL, 2*NL, 3*NL];
        rows_top_R = [1, NR + 1, 2*NR + 1];

        % External BC: stress-free upper surface of the first layer.
        L(ind_L(rows_top_L), :) = 0;
        L(ind_L(rows_top_L), ind_L) = S_layers{n}(rows_top_L, :);

        % Continuity BC: tractions between the first and second layers.
        L(ind_L(rows_bot_L), :) = 0;
        L(ind_L(rows_bot_L), ind_L) = S_layers{n}(rows_bot_L, :);
        L(ind_L(rows_bot_L), ind_R) = -S_layers{n + 1}(rows_top_R, :);

    elseif n == num_of_layers % Last layer
        NL = Ns(n - 1);
        NR = Ns(n);
        ind_L = get_layer_ind(Ns, n - 1);
        ind_R = get_layer_ind(Ns, n);
        rows_top_R = [1, NR + 1, 2*NR + 1];
        rows_bot_R = [NR, 2*NR, 3*NR];

        % Displacement BC: continuity between the previous and current layers.
        L(ind_R(rows_top_R), :) = 0;
        for comp = 1:3
            L(ind_R(rows_top_R(comp)), ind_L) = get_interface_BC(Is{n - 1}(NL, :), NL, comp);
            L(ind_R(rows_top_R(comp)), ind_R) = get_interface_BC(-Is{n}(1, :), NR, comp);
        end

        % External BC: stress-free lower surface of the last layer.
        L(ind_R(rows_bot_R), :) = 0;
        L(ind_R(rows_bot_R), ind_R) = S_layers{n}(rows_bot_R, :);

    else % Middle layers
        NL = Ns(n - 1);
        NM = Ns(n);
        NR = Ns(n + 1);
        ind_L = get_layer_ind(Ns, n - 1);
        ind_M = get_layer_ind(Ns, n);
        ind_R = get_layer_ind(Ns, n + 1);
        rows_top_M = [1, NM + 1, 2*NM + 1];
        rows_bot_M = [NM, 2*NM, 3*NM];
        rows_top_R = [1, NR + 1, 2*NR + 1];

        % Displacement BC: continuity between the previous and current layers.
        L(ind_M(rows_top_M), :) = 0;
        for comp = 1:3
            L(ind_M(rows_top_M(comp)), ind_L) = get_interface_BC(Is{n - 1}(NL, :), NL, comp);
            L(ind_M(rows_top_M(comp)), ind_M) = get_interface_BC(-Is{n}(1, :), NM, comp);
        end

        % Continuity BC: tractions between the current and next layers.
        L(ind_M(rows_bot_M), :) = 0;
        L(ind_M(rows_bot_M), ind_M) = S_layers{n}(rows_bot_M, :);
        L(ind_M(rows_bot_M), ind_R) = -S_layers{n + 1}(rows_top_R, :);
    end
end
end

function I_bc = get_interface_BC(I_n, N, Pos)
I_bc = zeros(1, 3*N);
I_bc(1, (N*(Pos-1) + 1):N*Pos) = I_n;
end

function ind = get_layer_ind(Ns, n)
ind = (3*sum(Ns(1:(n-1))) + 1):3*sum(Ns(1:n));
end

function [U11, U1n, U21, U2n, U31, U3n] = get_ML_border_displacements(U, Ns)
% Top surface of the first layer
U11 = U(1);
U21 = U(Ns(1) + 1);
U31 = U(2*Ns(1) + 1);

% Bottom surface of the last layer
n = length(Ns);
ind0 = 3*sum(Ns(1:(n-1)));
U1n = U(ind0 + Ns(n));
U2n = U(ind0 + 2*Ns(n));
U3n = U(ind0 + 3*Ns(n));
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
