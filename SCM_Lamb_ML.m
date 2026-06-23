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
% Thickness of plate 
h = 0.5e-3; % [m]

%% Material properties and general parameters
% Mechanical properties of layer
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
    M = get_Lamb_M_matrix(Ns(n), densities(n)); % It is possible to exclude M matrix calculation from cycle.
    M_layers{n} = set_Lamb_M_BC(M, Ns(n)); % Boundary conditions for M matrix
end
M = assemble_Lamb_ML_M_matrix(M_layers, Ns);

% Init output phase velocities
total_N = sum(Ns);
total_dof = 2*total_N;
freqs = NaN(length(wavenumber), total_dof);
PV = NaN(length(wavenumber), total_dof);
PV_S = NaN(length(wavenumber), total_dof);
PV_AS = NaN(length(wavenumber), total_dof);
k_plot =  NaN(length(wavenumber), total_dof); 

% Calculation dispersion curves
timestart = tic;

n_wavenumber = length(wavenumber);
f = waitbar(0, 'Calculations.', 'Name', 'Calculations...');

for i = 1:length(wavenumber)
    % calculating L and S matrices
    k = wavenumber(i);
    L_layers = cell(1, num_of_layers);
    S_layers = cell(1, num_of_layers);

    for n = 1:num_of_layers
        L_layers{n} = get_Lamb_L_matrix(c{n}, k, I{n}, DM1{n}, DM2{n});
        S_layers{n} = get_Lamb_S_matrix(c{n}, k, I{n}, DM1{n});
    end

    % Setting boundary and interface conditions
    % Boundary conditions for L matrix
    L = assemble_Lamb_ML_L_matrix(L_layers, S_layers, I, Ns);

    % Calculation eigenvalues and eigenvectors of equation L*U = w2*M*U
    % Balancing
    [T1, T2] = balance2(L, M);
    L_bal = T1*L*T2;
    M_bal = T1*M*T2;
    % Lamb mode decoupled
    [U, w2] = eig(L_bal, M_bal, 'qz');
    U = T2 * U;
    w2 = diag(w2);
    [w2, inds] = sort(w2);
    w = real(sqrt(w2));
    Vp = w./k;
    Vp = Vp';
    fs = w./(2*pi);
    freqs(i, :) = fs/1000;
    PV(i, :) = Vp;

    % Lamb mode separation
    % U = real(U);
    % U = U(:, inds);


    for j = 1:total_dof
        % The value we supposed to be zero or an infitisinal value for approximately equal values
        % This value is used instead of zero to avoid computational errors
        z_tol = 1e-9;
        % Modes separation
        if Vp(j) ~= 0 && isfinite(Vp(j)) % Skip values equal to zero or Inf
            % Values of normalized displacements
            [U11, U1n, U31, U3n] = get_Lamb_ML_border_displacements(U(:, j), Ns);
            
            if (abs(abs(U11 - U1n)) <= z_tol) && (abs(abs(U31) - abs(U3n)) <= z_tol) && (abs(U1n) >= z_tol || (abs(U31 + U3n) <= z_tol))
                % Symmetric mode:
                % Interface U1: displacements on the borders are equal;
                % Interface U3: modulus of displacements on the borbers are equal;
                % Displacements on the borders of U1 should be greater zero
                % or displacements on the borders of U3 should be opposite in sign
                PV_S(i, j) = Vp(j);
            else
                % Otherwise the mode is antisymmetric
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

fprintf('SCM calculation took %.3f s\n', toc(timestart));

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
title('All Lamb modes');



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

function L = get_Lamb_L_matrix(c, k, I, D, D2)
L11 = c(5,5)*D2 - c(1,1)*I*k^2 + 2*c(1,5)*k*1i*D;
L13 = c(3,5)*D2 - c(1,5)*I*k^2 + (c(1,3) + c(5,5))*k*1i*D;
L31 = L13;
L33 = c(3,3)*D2 - c(5,5)*I*k^2 + 2*c(3,5)*k*1i*D;
L = [L11 L13; ...
     L31 L33];
end

function S = get_Lamb_S_matrix(c, k, I, D)
S1 = c(1,3)*k*1i*I + c(3,5)*D;
S3 = c(3,5)*k*1i*I + c(3,3)*D;
S7 = c(1,5)*k*1i*I + c(5,5)*D;
S9 = c(5,5)*k*1i*I + c(3,5)*D;
S = [S1 S3; ...
     S7 S9];
end

function M = get_Lamb_M_matrix(N, p)
M = eye(N*2).*(-p);
end

function L = set_Lamb_L_BC (L, S, N)
L(1, 1:2*N) = S(1, 1:2*N);
L(N, 1:2*N) = S(N, 1:2*N);
L(N+1, 1:2*N) = S(N+1, 1:2*N);
L(2*N, 1:2*N) = S(2*N, 1:2*N);
end

function M = set_Lamb_M_BC(M, N)
M(1, 1) = 0;
M(N, N) = 0;
M((N+1), (N+1)) = 0;
M(2*N, 2*N) = 0;
end

function M = assemble_Lamb_ML_M_matrix(M_layers, Ns)
total_dof = 2*sum(Ns);
M = zeros(total_dof, total_dof);

for n = 1:length(Ns)
    ind = (2*sum(Ns(1:(n-1)))+1):2*sum(Ns(1:n));
    M(ind, ind) = M_layers{n};
end
end

function L = assemble_Lamb_ML_L_matrix(L_layers, S_layers, Is, Ns)
num_of_layers = length(Ns);
total_dof = 2*sum(Ns);
L = zeros(total_dof, total_dof);

if num_of_layers == 1
    L = set_Lamb_L_BC(L_layers{1}, S_layers{1}, Ns(1));
    return;
end

for n = 1:num_of_layers
    if n == 1 % First layer
        NL = Ns(n); NR = Ns(n + 1);
        % Lamb mode
        buff_L = L_layers{n};
        buff_R = zeros(2*NL, 2*NR);
        % External BC
        buff_L(1, :) = S_layers{n}(1, :);
        buff_L((NL), :) = S_layers{n}((NL+1), :);
        % Continuity BC
        % 1st layer
        buff_L((NL+1), :) = S_layers{n}((NL), :);
        buff_L((2*NL), :) = S_layers{n}((2*NL), :);
        % 2nd layer
        buff_R((NL+1), :) = -S_layers{(n+1)}((1), :);
        buff_R((2*NL), :) = -S_layers{(n+1)}((NR+1), :);
        % Apply first layer matrix
        L(1:2*sum(Ns(1:(n))), 1:2*sum(Ns(1:(n+1)))) = horzcat(buff_L, buff_R);

    elseif n == num_of_layers % Last layer
        NL = Ns(n - 1); NR = Ns(n);
        % Lamb mode
        buff_R = L_layers{n};
        buff_L = zeros(2*NR, 2*NL);
        % External BC
        buff_R((NR+1), :) = S_layers{n}((NR), :);
        buff_R((2*NR), :) = S_layers{n}((2*NR), :);
        % Displacement BC
        % (n-1)-th layer
        buff_L(1, :) = get_Lamb_interface_BC(Is{(n-1)}(NL, :), NL, 1);
        buff_L(NR, :) = get_Lamb_interface_BC(Is{(n-1)}(NL, :), NL, 2);
        % n-th layer
        buff_R(1, :) = get_Lamb_interface_BC(-Is{n}(1, :), NR, 1);
        buff_R(NR, :) = get_Lamb_interface_BC(-Is{n}(1, :), NR, 2);
        % Apply last layer matrix
        L((2*sum(Ns(1:(n - 1)))+1):2*sum(Ns(1:n)), (2*sum(Ns(1:(n - 2)))+1):(2*sum(Ns(1:n)))) = horzcat(buff_L, buff_R);

    else % Middle layers
        NL = Ns(n-1);
        NM = Ns(n);
        NR = Ns(n+1);
        % Lamb mode
        buff_L = zeros(2*NM, 2*NL);
        buff_M = L_layers{n};
        buff_R = zeros(2*NM, 2*NR);
        % Continuity BC
        % n-th layer
        buff_M((NM+1), :) = S_layers{n}(NM, :);
        buff_M((2*NM), :) = S_layers{n}((2*NM), :);
        % (n+1)-th layer
        buff_R((NM+1), :) = -S_layers{(n+1)}((1), :);
        buff_R((2*NM), :) = -S_layers{(n+1)}((NR+1), :);
        % Displacement BC
        % (n-1)-th layer
        buff_L(1, :) = get_Lamb_interface_BC(Is{(n-1)}((NL), :), NL, 1);
        buff_L((NM), :) = get_Lamb_interface_BC(Is{(n-1)}((NL), :), NL, 2);
        % n-th layer
        buff_M(1, :) = get_Lamb_interface_BC(-Is{n}(1, :), NM, 1);
        buff_M((NM), :) = get_Lamb_interface_BC(-Is{n}(1, :), NM, 2);
        % Apply n-th layer matrix
        L((2*sum(Ns(1:(n - 1)))+1):2*sum(Ns(1:n)), (2*sum(Ns(1:(n-2)))+1):2*sum(Ns(1:(n+1)))) = horzcat(buff_L, buff_M, buff_R);
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

