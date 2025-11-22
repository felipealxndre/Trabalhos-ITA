clear
close all
clc

global g
global aircraft

g = 9.80665;

% Dados do Dornier 328 (Relatório Pág. 6)
aircraft = struct('S', 40, 'c', 2.04, 'b', 20.8, ...
    'm', 10500, ...         % Massa base (Ajustar depois para Voo 1 e 2)
    'Ixx', 1.03e5, 'Iyy', 1.58e5, 'Izz', 2.40e5, 'Ixz', 1.25e4, ...
    'i_p_deg', -1.5, ...    % Incidência do motor
    'x_p', 1.50, 'z_p', -0.55, ... % Posição do motor
    'Tmax', 27350, ...      % Empuxo máx ao nível do mar
    'n_rho', 0.8, ...       % Expoente da densidade
    'n_V', -1, ...          % Expoente da velocidade (Turboélice)
    'V_ref', 120.49);       % Velocidade de referência

psidot_deg_s_eq = 0; % Voo reto e nivelado

% Condição de Cruzeiro (Relatório Pág. 6)
% NOTA: Para o trabalho, você precisará criar trim_par(2) e (3) com os dados
% de peso e altitude específicos do Voo 1 e Voo 2.
trim_par(1) = struct('V', 120.49, 'h', 5486.4, 'gamma_deg', 0, ...
    'thetadot_deg_s', 0, 'psidot_deg_s', psidot_deg_s_eq);

% Opções de otimização
options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);

%--------------------------------------------------------------------------
% Trimagem

% Inicializa output para U_eq com tamanho 4 (sem i_t)
trim_output(1) = struct('X_eq',zeros(12,1),'U_eq',zeros(4,1),'Y_eq',zeros(11,1));

for i_cond = 1:1
    % Vetor de chute inicial (x)
    % x = [V; alpha; q; theta; phi; p; r; psi; throttle; delta_e; delta_a; delta_r]
    x_eq_0 = zeros(12,1);
    x_eq_0(1) = trim_par(i_cond).V;
    x_eq_0(2) = 2;   % Chute inicial alpha (graus)
    x_eq_0(4) = 2;   % Chute inicial theta (graus)
    x_eq_0(9) = 0.4; % Chute inicial manete (0 a 1)
    
    % Função fsolve chama trim_function
    x_eq = fsolve(@trim_function, x_eq_0, options, trim_par, i_cond);
    
    [~,X_eq,U_eq,Y_eq] = trim_function(x_eq, trim_par, i_cond);
    
    trim_output(i_cond).X_eq = X_eq;
    trim_output(i_cond).U_eq = U_eq;
    trim_output(i_cond).Y_eq = Y_eq;

    fprintf('----- CONDIÇÃO DE VOO %d -----\n\n',i_cond);
    fprintf('   %-10s = %10.2f %-4s\n','V',X_eq(1),'m/s');
    fprintf('   %-10s = %10.4f %-4s\n','alpha',X_eq(2),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','theta',X_eq(4),'deg');
    fprintf('   %-10s = %10.1f %-4s\n','h',X_eq(5),'m');
    fprintf('\n');
    fprintf('   %-10s = %10.4f %-4s\n','beta',X_eq(7),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','phi',X_eq(8),'deg');
    fprintf('\n');
    fprintf('   %-10s = %10.2f %-4s\n','throttle',U_eq(1)*100,'%');
    fprintf('   %-10s = %10.4f %-4s\n','delta_e',U_eq(2),'deg'); % Índice mudou para 2
    fprintf('   %-10s = %10.4f %-4s\n','delta_a',U_eq(3),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','delta_r',U_eq(4),'deg');
    fprintf('\n');
end

save trim_output.mat trim_output


%% --- SIMULAÇÃO DO PERÍODO CURTO (Doublet de Profundor) ---

% 1. Configuração Inicial
i_cond = 1;
X_trim = trim_output(i_cond).X_eq;
U_trim = trim_output(i_cond).U_eq;

% 2. Definição do Tempo de Simulação
t_final = 10; % O período curto é rápido, 10s é suficiente
dt = 0.01;
t_span = 0:dt:t_final;

% 3. Execução da Simulação
% Usamos uma função anônima para injetar o doublet no controle U(2) = delta_e
% O doublet ocorre entre t=1s e t=2s
[T_sim, X_sim] = ode45(@(t,x) dynamics_with_input(t, x, U_trim, i_cond), t_span, X_trim);

% 4. Plotagem dos Resultados
figure('Name', 'Resposta do Período Curto');

subplot(3,1,1);
plot(T_sim, X_sim(:,2)); % Alpha
grid on; ylabel('\alpha [deg]'); title('Ângulo de Ataque');

subplot(3,1,2);
plot(T_sim, X_sim(:,3)); % q
grid on; ylabel('q [deg/s]'); title('Taxa de Arfagem');

subplot(3,1,3);
plot(T_sim, X_sim(:,1)); % V
grid on; ylabel('V [m/s]'); xlabel('Tempo [s]'); title('Velocidade');

%% --- FUNÇÃO AUXILIAR PARA O INPUT (Coloque no final do script) ---
function [dx] = dynamics_with_input(t, x, U_trim, i_cond)
    % Copia o controle trimado
    U_current = U_trim;
    
    amplitude = 2; % Amplitude do pulso em graus
    t_start = 1.0;
    duration = 1.0; % Duração total do doublet
    
    if t >= t_start && t < (t_start + duration/2)
        U_current(2) = U_trim(2) - amplitude; % Pulso negativo (cabrar)
    elseif t >= (t_start + duration/2) && t < (t_start + duration)
        U_current(2) = U_trim(2) + amplitude; % Pulso positivo (picar)
    end
    
    % Chama a dinâmica original
    [dx, ~] = dynamics(t, x, U_current, i_cond);
end

%% --- ANÁLISE LINEAR (Autovalores) ---

i_cond = 1;
X_eq = trim_output(i_cond).X_eq;
U_eq = trim_output(i_cond).U_eq;

nX = length(X_eq);
perturbacao = 1e-5;
A = zeros(nX, nX);

% Cálculo Numérico da Matriz Jacobiana A
for k = 1:nX
    X_plus = X_eq; X_plus(k) = X_plus(k) + perturbacao;
    X_minus = X_eq; X_minus(k) = X_minus(k) - perturbacao;
    
    [dX_plus, ~] = dynamics(0, X_plus, U_eq, i_cond);
    [dX_minus, ~] = dynamics(0, X_minus, U_eq, i_cond);
    
    A(:, k) = (dX_plus - dX_minus) / (2 * perturbacao);
end

% Separar modos Longitudinais
% Estados Longitudinais Típicos: V(1), alpha(2), q(3), theta(4)
idx_long = [1, 2, 3, 4]; 
A_long = A(idx_long, idx_long);

% Calcular Autovalores
autovalores = eig(A_long);

fprintf('\n--- AUTOVALORES LONGITUDINAIS ---\n');
disp(autovalores);

% Função 'damp' do Matlab já calcula wn e zeta para você
fprintf('\n--- CARACTERÍSTICAS DOS MODOS ---\n');
damp(autovalores);
