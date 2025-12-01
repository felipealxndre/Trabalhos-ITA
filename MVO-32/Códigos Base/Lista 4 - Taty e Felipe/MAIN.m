clear
close all
clc

global g
global aircraft
create_aircraft;

g = 9.80665;

trim_par = struct('V',228.8138886191133,'h',11000);

options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);

%% Ex.2 (Lista 3) – voo reto nivelado em cruzeiro
% Este equilíbrio será o ponto de operação base da Lista 4.

x_eq_0_Ex2 = zeros(5,1);
x_eq_Ex2 = fsolve(@trim_function_Ex2,x_eq_0_Ex2,options,trim_par);
[~,X_eq_Ex2,U_eq_Ex2,Y_eq_Ex2] = trim_function_Ex2(x_eq_Ex2,trim_par);

fprintf('\n----- EX.2 -----\n\n');
fprintf('   %-12s = %8.2f %s\n','V',X_eq_Ex2(1),'m/s');
fprintf('   %-12s = %8.3f %s\n','alpha',X_eq_Ex2(2),'deg');
fprintf('   %-12s = %8.3f %s\n','q',X_eq_Ex2(3),'deg/s');
fprintf('   %-12s = %8.3f %s\n','theta',X_eq_Ex2(4),'deg');
fprintf('   %-12s = %8.1f %s\n','h',X_eq_Ex2(5),'m');
fprintf('   %-12s = %8.3f %s\n','beta',X_eq_Ex2(7),'deg');
fprintf('   %-12s = %8.3f %s\n','phi',X_eq_Ex2(8),'deg');
fprintf('   %-12s = %8.3f %s\n','p',X_eq_Ex2(9),'deg/s');
fprintf('   %-12s = %8.3f %s\n','r',X_eq_Ex2(10),'deg/s');
fprintf('\n');
fprintf('   %-12s = %8.2f %s\n','delta1',U_eq_Ex2(1),'');
fprintf('   %-12s = %8.2f %s\n','delta2',U_eq_Ex2(2),'');
fprintf('   %-12s = %8.3f %s\n','i_t',U_eq_Ex2(3),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_e',U_eq_Ex2(4),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_a',U_eq_Ex2(5),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_r',U_eq_Ex2(6),'deg');
fprintf('\n');
fprintf('   %-12s = %8.3f %s\n','gamma',Y_eq_Ex2(1),'deg');
fprintf('   %-12s = %8.1f %s\n','T1',Y_eq_Ex2(2),'N');
fprintf('   %-12s = %8.1f %s\n','T2',Y_eq_Ex2(3),'N');
fprintf('   %-12s = %8.3f %s\n','Mach',Y_eq_Ex2(4),'');
fprintf('   %-12s = %8.4f %s\n','C_D',Y_eq_Ex2(5),'');
fprintf('   %-12s = %8.4f %s\n','C_L',Y_eq_Ex2(6),'');
fprintf('   %-12s = %8.4f %s\n','C_m',Y_eq_Ex2(7),'');
fprintf('   %-12s = %8.4f %s\n','C_Y',Y_eq_Ex2(8),'');
fprintf('   %-12s = %8.4f %s\n','C_l',Y_eq_Ex2(9),'');
fprintf('   %-12s = %8.4f %s\n','C_n',Y_eq_Ex2(10),'');
fprintf('   %-12s = %8.4f %s\n','rho',Y_eq_Ex2(11),'kg/m^3');
fprintf('   %-12s = %8.1f %s\n','qbar',Y_eq_Ex2(12),'N/m^2');
fprintf('\n');

%% Ex.3 (Lista 3) – linearização em torno do equilíbrio de Ex.2
% ESTE É O MODELO QUE A LISTA 4 PEDE ("modelo de aeronave desenvolvido na Lista 3")

nX_Ex3 = length(X_eq_Ex2);
nU_Ex3 = length(U_eq_Ex2);
nY_Ex3 = length(Y_eq_Ex2);

lin_output_Ex3 = struct('A',zeros(nX_Ex3,nX_Ex3),'B',zeros(nX_Ex3,nU_Ex3), ...
                        'C',zeros(nY_Ex3,nX_Ex3),'D',zeros(nY_Ex3,nU_Ex3));

delta_val_Ex3 = 1e-5;

A_Ex3 = zeros(nX_Ex3,nX_Ex3);
C_Ex3 = zeros(nY_Ex3,nX_Ex3);
for j=1:nX_Ex3
    dX = zeros(nX_Ex3,1);
    dX(j) = delta_val_Ex3;
    [Xdot_plus, Y_plus]= dynamics(0, X_eq_Ex2 + dX, U_eq_Ex2);
    [Xdot_minus, Y_minus]= dynamics(0, X_eq_Ex2 - dX, U_eq_Ex2);
    A_Ex3(:, j) = (Xdot_plus - Xdot_minus)/(2*delta_val_Ex3);
    C_Ex3(:, j) = (Y_plus - Y_minus)/(2*delta_val_Ex3);
end

B_Ex3 = zeros(nX_Ex3,nU_Ex3);
D_Ex3 = zeros(nY_Ex3,nU_Ex3);
for j=1:nU_Ex3
    dU = zeros(nU_Ex3,1);
    dU(j) = delta_val_Ex3;
    [Xdot_plus, Y_plus]= dynamics(0, X_eq_Ex2, U_eq_Ex2 + dU);
    [Xdot_minus, Y_minus]= dynamics(0, X_eq_Ex2, U_eq_Ex2 - dU);
    B_Ex3(:, j) = (Xdot_plus - Xdot_minus)/(2*delta_val_Ex3);
    D_Ex3(:, j) = (Y_plus - Y_minus)/(2*delta_val_Ex3);
end

lin_output_Ex3.A = A_Ex3;
lin_output_Ex3.B = B_Ex3;
lin_output_Ex3.C = C_Ex3;
lin_output_Ex3.D = D_Ex3;

sel_states_Ex3 = 1:12;
[eigvec_Ex3,eigval_Ex3] = eig(lin_output_Ex3.A(sel_states_Ex3,sel_states_Ex3));
eigval_Ex3 = diag(eigval_Ex3);
[eigval_Ex3,i_sort_Ex3] = sort(eigval_Ex3);
eigvec_Ex3 = eigvec_Ex3(:,i_sort_Ex3);

fprintf('\n----- EX.3 (Lista 3) - autovalores do modelo base -----\n\n');
damp(eigval_Ex3)

%% L4 – Modelo base para todos os exercícios (salvar para usar em scripts da Lista 4)
% Vamos guardar o modelo linear completo em malha aberta da Lista 3
% para usar nos scripts de:
% - L4 Ex.1 (projeto SAS longitudinal por lugar das raízes)
% - L4 Ex.2 (simulação do doublet com SAS)
% - L4 Ex.3 (modelo sem dinâmica de motores)

lin_output_L4 = lin_output_Ex3;   % só renomeando para ficar semântico
X_eq_L4      = X_eq_Ex2;
U_eq_L4      = U_eq_Ex2;
Y_eq_L4      = Y_eq_Ex2;

save('lin_L4_modelo.mat','lin_output_L4','X_eq_L4','U_eq_L4','Y_eq_L4');

fprintf('\nArquivo lin_L4_modelo.mat salvo com o modelo linear completo em malha aberta.\n');
fprintf('Use este arquivo nos scripts da Lista 4 (projeto do SAS, simulações etc.).\n\n');

%% L4 – Ex.3 – inclusão das dinâmicas de primeira ordem dos motores
% Equação:  d(delta_i)/dt = -(1/tau_throttle)*delta_i + (1/tau_throttle)*u_i
% i = 1 (motor esquerdo), i = 2 (motor direito)
% Aqui partimos de xdot = A x + B u, y = C x + D u,
% com u = [delta1; delta2; i_t; delta_e; delta_a; delta_r].
%
% Novo modelo aumentado:
% x_aug = [x; delta1; delta2]
% u_aug = [u1; u2; i_t; delta_e; delta_a; delta_r] (mesma ordem)
%
% x_dot   = A x + B_delta * [delta1; delta2] + B_rest * u_rest
% delta1_dot = -(1/tau)*delta1 + (1/tau)*u1
% delta2_dot = -(1/tau)*delta2 + (1/tau)*u2

tau_throttle = 0.1;  % s

A_base = lin_output_L4.A;
B_base = lin_output_L4.B;
C_base = lin_output_L4.C;
D_base = lin_output_L4.D;

[nx,nu] = size(B_base);
ny = size(C_base,1);

% Separar colunas de B associadas aos motores e aos demais controles
% u = [delta1; delta2; i_t; delta_e; delta_a; delta_r]
B_delta = B_base(:,1:2);      % efeito de delta1, delta2
B_rest  = B_base(:,3:end);    % efeito de i_t, delta_e, delta_a, delta_r

% Matriz A aumentada (nx+2 estados)
A_throttle = zeros(nx+2,nx+2);
A_throttle(1:nx,1:nx) = A_base;
A_throttle(1:nx,nx+1:nx+2) = B_delta;        % x_dot depende de delta1, delta2
A_throttle(nx+1,nx+1) = -1/tau_throttle;     % delta1_dot
A_throttle(nx+2,nx+2) = -1/tau_throttle;     % delta2_dot

% Matriz B aumentada (mesmo número de entradas)
B_throttle = zeros(nx+2,nu);
% primeira linha até nx: x_dot só depende de i_t, delta_e, delta_a, delta_r
B_throttle(1:nx,3:nu) = B_rest;
% estados delta1, delta2 com dinâmicas de primeira ordem
B_throttle(nx+1,1) = 1/tau_throttle;   % u1 -> delta1_dot
B_throttle(nx+2,2) = 1/tau_throttle;   % u2 -> delta2_dot

% Matriz C aumentada: saídas originais só enxergam x (por simplicidade)
C_throttle = [C_base, zeros(ny,2)];

% Matriz D aumentada: mantemos a mesma
D_throttle = D_base;

% Autovalores sem dinâmica de motores
lambda_sem_throttle = eig(A_base);

% Autovalores com dinâmica de motores
lambda_com_throttle = eig(A_throttle);

fprintf('\n----- L4 Ex.3 - autovalores sem dinâmica dos motores -----\n\n');
damp(lambda_sem_throttle)

fprintf('\n----- L4 Ex.3 - autovalores COM dinâmica dos motores -----\n\n');
damp(lambda_com_throttle)

% Opcional: salvar o modelo aumentado para uso em outros scripts
lin_output_L4_throttle = struct('A',A_throttle,'B',B_throttle, ...
                                'C',C_throttle,'D',D_throttle);

save('lin_L4_modelo_throttle.mat','lin_output_L4_throttle','X_eq_L4','U_eq_L4','Y_eq_L4');

fprintf('\nArquivo lin_L4_modelo_throttle.mat salvo com o modelo com dinâmicas de motores.\n');
fprintf('Note que aparecem dois autovalores reais aproximadamente iguais a -1/tau_throttle = -10.\n\n');

