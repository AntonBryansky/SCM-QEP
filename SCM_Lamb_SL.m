%% Startup cleanup
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
freq_limit = 1e6;           % Hz
% Wavenumber limit 
WN_limit = 4000.0;          % rad/m for plotting y-axis
WN_amount = 500;            % Amount of wavenumber scanning points
wavenumber = 0:(WN_limit/WN_amount):WN_limit;
wavenumber(1) = 1;          % avoid zero values
PV_limit = 12000.0;         % Limit of the phase velocity for plots, used for plot


%% Transform stiffness matrix relatively to the propagation direction (angle)
% Angle between the main direction and the plate orientations
beta = orientation - propagation_angle;
% Transform the stiffness matrix to the propagation direction
c = transform_stiffness_matrix(C_CFRP, beta);

%% Calculation
% Chebyshev differentiation matrices
% Calculation of the differentiation matrices for domain [-1 1]
[s, DMs] = chebdif(N, 2);
% Scaling coefficient for domain [-h/2 h/2]
scale = 2.0 / h;
DM1 = DMs(:,:,1) * scale;
DM2 = DMs(:,:,2) * scale^2;

% Physical thickness coordinate, m
z = (h/2) * s; 
z_mm = z * 1e3;             % mm
I = eye(N);                 % identity matrix

% Independent matrices
M = get_Lamb_M_matrix(N, density);      % M matrix does not depend on the wavenumber
M = set_Lamb_M_BC(M, N);                % Setting boundary conditions for M matrix

% Init output data storage
freqs_output = NaN(length(wavenumber), N*2);
PV_output = NaN(length(wavenumber),N*2);
PV_S_output = NaN(length(wavenumber),N*2);
PV_AS_output = NaN(length(wavenumber),N*2);
k_output =  NaN(length(wavenumber),N*2);

% Calculation dispersion curves
timestart = tic; % intitialize timer
f = waitbar(0, 'Calculations.', 'Name', 'Calculations...');

n_wavenumber = length(wavenumber);

for i = 1:n_wavenumber
    k = wavenumber(i);
    % calculating L and S matrices
    L = get_Lamb_L_matrix(c, k, I, DM1, DM2);
    S = get_Lamb_S_matrix(c, k, I, DM1);
    
    % Setting boundary and interface conditions
    % Boundary conditions for L matrix
    L = set_Lamb_L_BC(L, S, N);
    
    % Calculation eigenvalues and eigenvectors of equation L*U = w2*M*U
    % Lamb mode is decoupled
    [U, w2] = eig(L, M, 'qz');
    w2 = diag(w2);              % squared angular frequency
    [w2, inds] = sort(w2);
    w = real(sqrt(w2));         % angular frequency, 1/s
    Vp = w./k;                  % phase velocity, m/s
    Vp = Vp';
    freqs_output(i, :) = w./(2*pi); % frequency, Hz
    PV_output(i, :) = Vp;
    % Lamb mode separation
    % Displacement profile components
    U = real(U);
    U = U(:, inds);
    U1 = U(1:N, :);
    U3 = U((N+1):2*N, :);
    for j = 1:2*N
        % The value we supposed to be zero or an infitisinal value for approximately equal values
        % This value is used instead of zero to avoid computational errors
        err_coeff = 1e-7;
        % Modes separation
        if Vp(j) ~= 0 && isfinite(Vp(j)) % Skip values equal to zero or Inf
            % Values of normalized displacements
            U11 = U1(1, j); U1n = U1(end, j);
            U31 = U3(1, j); U3n = U3(end, j);
            if (abs(abs(U11 - U1n)) <= err_coeff) && (abs(abs(U31) - abs(U3n)) <= err_coeff) && (abs(U1n) >= err_coeff || (abs(U31 + U3n) <= err_coeff))
                % Symmetric mode:
                % Interface U1: displacements on the borders are equal;
                % Interface U3: modulus of displacements on the borbers are equal;
                % Displacements on the borders of U1 should be greater zero
                % or displacements on the borders of U3 should be opposite in sign
                PV_S_output(i, j) = Vp(j);
            else
                % Otherwise the mode is antisymmetric
                PV_AS_output(i, j) = Vp(j);
            end
        end
    end
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
waitbar(100, f, 'Calculations completed!');
delete(f);

fprintf('SCM calculation took %.3f s\n', toc(timestart));

%% Plot
% Phase velocity versus frequency plot

figure('Color','w');
p1 = plot(freqs_output*1e-3, PV_S_output, '.r');
p1 = p1(1);
ylim([0 PV_limit]);
xlim([0 1000]);
hold on;
p2 = plot(freqs_output*1e-3, PV_AS_output, '.b');
p2 = p2(1);
ylim([0 PV_limit]);
xlim([0 freq_limit*1e-3]);
xlabel('Frequency, kHz');
ylabel('Phase velocity, m/s');
legend([p1, p2], {'SM', 'ASM'}, 'Location', 'best');



% Wavenumenber versus frequency plot

Kmat = repmat(wavenumber(:), 1, 2*N);   % each row has the same input k
% Use your already separated symmetric/antisymmetric masks
K_S  = Kmat;
K_AS = Kmat;
K_S(~isfinite(PV_S_output))   = NaN;
K_AS(~isfinite(PV_AS_output)) = NaN;

figure('Color','w');
p1 = plot(freqs_output*1e-3, K_S,  '.r', 'MarkerSize', 4);
hold on;
p2 = plot(freqs_output*1e-3, K_AS, '.b', 'MarkerSize', 4);
xlim([0, freq_limit*1e-3]);
ylim([0, WN_limit]);
xlabel('Frequency (kHz)');
ylabel('Wavenumber Re(k_x) (1/m)');
legend([p1(1), p2(1)], {'SM', 'ASM'}, 'Location', 'best');
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
% if  c(1,6) == 0
%     c(1,6) = 1;
%     c(2,6) = 1;
%     c(3,6) = 1;
%     c(4,5) = 1;
% end
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
