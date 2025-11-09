clear
close all
clc

global g
global aircraft
aircraft = create_aircraft();

g = 9.80665;

trim_par = struct('V',228.8138886191133,'h',11000,'gamma_deg',0,...
    'thetadot_deg_s',0);

% https://www.mathworks.com/help/optim/ug/tolerances-and-stopping-criteria.html
options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);
% % In newer MATLAB versions, the following command might be necessary:
% options = optimoptions(@fsolve,'Display','iter','StepTolerance',1e-10,'FunctionTolerance',1e-10);
% % or:
% options = optimoptions(@fsolve,'Display','iter','TolX',1e-10,'TolFun',1e-10);

%% Ex.2

x_eq_0 = zeros(6,1);   
x_eq = fsolve(@trim_function,x_eq_0,options,trim_par);
[~,X_eq,U_eq,Y_eq] = trim_function(x_eq,trim_par);

fprintf('\n----- EX.2 -----\n\n');
fprintf('   %-12s = %8.2f %s\n','V',X_eq(1),'m/s');
fprintf('   %-12s = %8.3f %s\n','alpha',X_eq(2),'deg');
fprintf('   %-12s = %8.3f %s\n','q',X_eq(3),'deg/s');
fprintf('   %-12s = %8.3f %s\n','theta',X_eq(4),'deg');
fprintf('   %-12s = %8.1f %s\n','h',X_eq(5),'m');
fprintf('   %-12s = %8.3f %s\n','beta',X_eq(7),'deg');
fprintf('   %-12s = %8.3f %s\n','phi',X_eq(8),'deg');
fprintf('   %-12s = %8.3f %s\n','p',X_eq(9),'deg/s');
fprintf('   %-12s = %8.3f %s\n','r',X_eq(10),'deg/s');
fprintf('\n');
fprintf('   %-12s = %8.2f %s\n','delta1',U_eq(1),'%');
fprintf('   %-12s = %8.2f %s\n','delta2',U_eq(2),'%');
fprintf('   %-12s = %8.3f %s\n','i_t',U_eq(3),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_e',U_eq(4),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_a',U_eq(5),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_r',U_eq(6),'deg');
fprintf('\n');
fprintf('   %-12s = %8.3f %s\n','gamma',Y_eq(1),'deg');
fprintf('   %-12s = %8.1f %s\n','T1',Y_eq(2),'N');
fprintf('   %-12s = %8.1f %s\n','T2',Y_eq(3),'N');
fprintf('   %-12s = %8.3f %s\n','Mach',Y_eq(4),'');
fprintf('   %-12s = %8.4f %s\n','C_D',Y_eq(5),'');
fprintf('   %-12s = %8.4f %s\n','C_L',Y_eq(6),'');
fprintf('   %-12s = %8.4f %s\n','C_m',Y_eq(7),'');
fprintf('   %-12s = %8.4f %s\n','C_Y',Y_eq(8),'');
fprintf('   %-12s = %8.4f %s\n','C_l',Y_eq(9),'');
fprintf('   %-12s = %8.4f %s\n','C_n',Y_eq(10),'');
fprintf('   %-12s = %8.4f %s\n','rho',Y_eq(11),'kg/m^3');
fprintf('   %-12s = %8.1f %s\n','qbar',Y_eq(12),'N/m^2');
fprintf('\n');

%% Ex.3

% 1) Obter ponto de equilíbrio (usa dynamics 6-DoF)
%    trim_par: struct com V, h, gamma_deg, thetadot_deg_s
x_eq_0 = zeros(12,1);                 % [V,alpha,q,theta,h,x,beta,phi,p,r,psi,y] (deg/deg/s conforme seu modelo)
x_eq_0(1) = trim_par.V;               % chute inicial V
[X_res, ~, ~, ~] = fsolve(@trim_function, x_eq_0, options, trim_par);
[~, X_eq, U_eq, Y_eq] = trim_function(X_res, trim_par);

% 2) Dimensões
nX = length(X_eq);
nU = length(U_eq);
[~, Y_tmp] = dynamics(0, X_eq, U_eq);
nY = length(Y_tmp);

% 3) Deltas por unidade (evita misturar m, m/s, deg, deg/s e frações)
%    Ajuste mínimos conforme sua escala
delta_X = zeros(nX,1);
% estados: [V, alpha, q, theta, h, x, beta, phi, p, r, psi, y]
delta_X(1)  = max(1e-6, abs(X_eq(1))*1e-6);   % V [m/s]
delta_X(2)  = 1e-4;                            % alpha [deg]
delta_X(3)  = 1e-5;                            % q [deg/s]
delta_X(4)  = 1e-4;                            % theta [deg]
delta_X(5)  = max(1e-6, abs(X_eq(5))*1e-6);   % h [m]
delta_X(6)  = 1e-3;                            % x [m]
delta_X(7)  = 1e-4;                            % beta [deg]
delta_X(8)  = 1e-4;                            % phi [deg]
delta_X(9)  = 1e-5;                            % p [deg/s]
delta_X(10) = 1e-5;                            % r [deg/s]
delta_X(11) = 1e-4;                            % psi [deg]
delta_X(12) = 1e-3;                            % y [m]

% entradas U: [delta1, delta2, i_t, delta_e, delta_a, delta_r] (exemplo)
delta_U = zeros(nU,1);
delta_U(1) = 1e-4;   % throttle 1 (fração 0–1)
delta_U(2) = 1e-4;   % throttle 2 (fração 0–1)
delta_U(3) = 1e-4;   % i_t [deg]
delta_U(4) = 1e-4;   % delta_e [deg]
delta_U(5) = 1e-4;   % delta_a [deg]
delta_U(6) = 1e-4;   % delta_r [deg]

% 4) Matrizes A, B, C, D por diferenças centrais
A = zeros(nX,nX);  B = zeros(nX,nU);
C = zeros(nY,nX);  D = zeros(nY,nU);

% A e C (variação em X)
for j = 1:nX
    dX = zeros(nX,1); dX(j) = delta_X(j);
    [Xdot_p, Y_p] = dynamics(0, X_eq + dX, U_eq);
    [Xdot_m, Y_m] = dynamics(0, X_eq - dX, U_eq);
    A(:,j) = (Xdot_p - Xdot_m) / (2*delta_X(j));
    C(:,j) = (Y_p     - Y_m)    / (2*delta_X(j));
end

% B e D (variação em U)
for j = 1:nU
    dU = zeros(nU,1); dU(j) = delta_U(j);
    [Xdot_p, Y_p] = dynamics(0, X_eq, U_eq + dU);
    [Xdot_m, Y_m] = dynamics(0, X_eq, U_eq - dU);
    B(:,j) = (Xdot_p - Xdot_m) / (2*delta_U(j));
    D(:,j) = (Y_p     - Y_m)    / (2*delta_U(j));
end

lin_output = struct('A',A,'B',B,'C',C,'D',D, ...
                    'X_eq',X_eq,'U_eq',U_eq,'Y_eq',Y_eq, ...
                    'delta_X',delta_X,'delta_U',delta_U);
save lin_output.mat lin_output

% 5) Autovalores e amortecimento/frequência natural
[eigvec,eigmat] = eig(A);
eigval = diag(eigmat);
[~, idx] = sort(real(eigval), 'ascend');
eigval = eigval(idx);
eigvec = eigvec(:,idx);

fprintf('\n----- EX.3 -----\n\n');
damp(eigval)   % zeta e wn diretamente dos polos


%% Ex.4

x0 = [ X_eq(2); U_eq(1); U_eq(3); 0.0; 0.5; -1.0 ]; % chute: [alpha, δ1, i_t, δe, δa, δr]
[x_sol, ~] = fsolve(@trim_function_failure, x0, options, trim_par);
[~, X_eq_fail, U_eq_fail, Y_eq_fail] = trim_function_failure(x_sol, trim_par);

fprintf('\n----- EX.4 -----\n\n');
fprintf('   %-12s = %8.2f %s\n','V',X_eq_fail(1),'m/s');
fprintf('   %-12s = %8.3f %s\n','alpha',X_eq_fail(2),'deg');
fprintf('   %-12s = %8.3f %s\n','q',X_eq_fail(3),'deg/s');
fprintf('   %-12s = %8.3f %s\n','theta',X_eq_fail(4),'deg');
fprintf('   %-12s = %8.1f %s\n','h',X_eq_fail(5),'m');
fprintf('   %-12s = %8.3f %s\n','beta',X_eq_fail(7),'deg');
fprintf('   %-12s = %8.3f %s\n','phi',X_eq_fail(8),'deg');
fprintf('   %-12s = %8.3f %s\n','p',X_eq_fail(9),'deg/s');
fprintf('   %-12s = %8.3f %s\n','r',X_eq_fail(10),'deg/s');
fprintf('\n');
fprintf('   %-12s = %8.2f %s\n','delta1',U_eq_fail(1),'%');
fprintf('   %-12s = %8.2f %s\n','delta2',U_eq_fail(2),'%');
fprintf('   %-12s = %8.3f %s\n','i_t',U_eq_fail(3),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_e',U_eq_fail(4),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_a',U_eq_fail(5),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_r',U_eq_fail(6),'deg');
fprintf('\n');
fprintf('   %-12s = %8.3f %s\n','gamma',Y_eq_fail(1),'deg');
fprintf('   %-12s = %8.1f %s\n','T1',Y_eq_fail(2),'N');
fprintf('   %-12s = %8.1f %s\n','T2',Y_eq_fail(3),'N');
fprintf('   %-12s = %8.3f %s\n','Mach',Y_eq_fail(4),'');
fprintf('   %-12s = %8.4f %s\n','C_D',Y_eq_fail(5),'');
fprintf('   %-12s = %8.4f %s\n','C_L',Y_eq_fail(6),'');
fprintf('   %-12s = %8.4f %s\n','C_m',Y_eq_fail(7),'');
fprintf('   %-12s = %8.4f %s\n','C_Y',Y_eq_fail(8),'');
fprintf('   %-12s = %8.4f %s\n','C_l',Y_eq_fail(9),'');
fprintf('   %-12s = %8.4f %s\n','C_n',Y_eq_fail(10),'');
fprintf('   %-12s = %8.4f %s\n','rho',Y_eq_fail(11),'kg/m^3');
fprintf('   %-12s = %8.1f %s\n','qbar',Y_eq_fail(12),'N/m^2');
fprintf('\n');

%% Ex.5

%% Ex.6

%% Ex.7