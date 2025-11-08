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
x_eq_0(1) = trim_par.V;    
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

fprintf('\n----- EX.3 -----\n\n');

% Pega o Xdot de equilíbrio (deve ser [0; 0; 0; ...])
[Xdot_eq, ~] = dynamics(0, X_eq, U_eq);

% Inicializa a matriz A (12x12)
n_states = length(X_eq);
A = zeros(n_states, n_states);

% Define o tamanho da perturbação
perturb = 1e-6; 

% Loop para calcular cada coluna da Matriz A
for j = 1:n_states
    % Cria um vetor de perturbação
    X_pert = X_eq;
    
    % Perturba o j-ésimo estado
    X_pert(j) = X_pert(j) + perturb;
    
    % Calcula o Xdot com o estado perturbado
    [Xdot_pert, ~] = dynamics(0, X_pert, U_eq);
    
    % Calcula a derivada (coluna da matriz A) por diferenças finitas
    A_col = (Xdot_pert - Xdot_eq) / perturb;
    
    % Armazena a coluna na matriz A
    A(:, j) = A_col;
end

% Calcula autovalores (D_eig) e autovetores (V_eig)
[V_eig, D_eig] = eig(A);

% extrai os autovalores da diagonal de D_eig
eigvals = diag(D_eig);

fprintf('--- Análise dos Modos Naturais ---\n');
fprintf('%-20s | %-22s | %-12s | %-12s | %-10s\n', ...
    'Modo (Identifique!)', 'Autovalor (lambda)', 'Amort. (zeta)', 'Freq. (wn)', 'T. Const (tau)');
fprintf([repmat('-', 1, 84) '\n']);

% Tolerância para identificar autovalores nulos ou pares
tol = 1e-4; 
processed_mask = false(n_states, 1); % Para marcar autovalores já processados

for i = 1:n_states
    if processed_mask(i)
        continue; % Já foi processado como parte de um par
    end
    
    lambda = eigvals(i);
    
    if abs(lambda) < tol
        fprintf('%-20s | %-22s | %-12s | %-12s | %-10s\n', ...
            'Modo Ignorável (x,y,psi)', sprintf('%+.4f', real(lambda)), 'N/A', 'N/A', 'Infinito');
        processed_mask(i) = true;
        
    % Verifica se é um MODO REAL (Não-Oscilatório)
    elseif abs(imag(lambda)) < tol
        lambda = real(lambda); % Trata como puramente real
        tau = -1 / lambda;
        
        % Determina estabilidade
        if lambda > 0
            status = '(Instável)';
        else
            status = '(Estável)';
        end
        
        fprintf('%-20s | %-22s | %-12s | %-12s | %-10s\n', ...
            'Modo Real', sprintf('%+.4f %s', lambda, status), 'N/A', 'N/A', sprintf('%.3f s', tau));
        processed_mask(i) = true;

    % Verifica se é um MODO COMPLEXO (Oscilatório)
    else
        % Tenta encontrar o par conjugado
        conjugate_found = false;
        for k = i+1:n_states
            if abs(lambda - conj(eigvals(k))) < tol
                processed_mask(k) = true; % Marca o par como processado
                conjugate_found = true;
                break;
            end
        end
        
        % Se encontrou o par (ou se for o primeiro do par), calcula métricas
        if conjugate_found || ~processed_mask(i)
            sigma = real(lambda);
            omega_d = imag(lambda);
            omega_n = abs(lambda);
            zeta = -sigma / omega_n;
            
            % Determina estabilidade
            if sigma > 0
                status = '(Instável)';
            else
                status = '(Estável)';
            end

            fprintf('%-20s | %-22s | %-12s | %-12s | %-10s\n', ...
                'Modo Oscilatório', ...
                sprintf('%+.4f +/- %.4fj', sigma, abs(omega_d)), ...
                sprintf('%.4f %s', zeta, status), ...
                sprintf('%.4f rad/s', omega_n), ...
                'N/A');
            
            processed_mask(i) = true;
        end
    end
end
fprintf('\n');


%% Ex.4

%% Ex.5

%% Ex.6

%% Ex.7