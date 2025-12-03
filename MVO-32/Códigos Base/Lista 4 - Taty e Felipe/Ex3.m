clear
clc

% Gera modelo da Lista 3
[A,B,C,D,X_eq,U_eq,Y_eq] = gera_modelo_L3();

fprintf('\n========== L4 - EXERCÍCIO 3: DINÂMICAS DOS MOTORES ==========\n');

% Parâmetro do atuador de throttle
tau_throttle = 0.1;    % s

A_base = A;
B_base = B;
C_base = C;
D_base = D;

[nx,nu] = size(B_base);
ny = size(C_base,1);

% Entradas: u = [delta1; delta2; i_t; delta_e; delta_a; delta_r]
B_delta = B_base(:,1:2);      % efeito de delta1, delta2
B_rest  = B_base(:,3:end);    % efeito de i_t, delta_e, delta_a, delta_r

% Matriz A aumentada (inclui estados delta1 e delta2)
A_throttle = zeros(nx+2,nx+2);
A_throttle(1:nx,1:nx) = A_base;
A_throttle(1:nx,nx+1:nx+2) = B_delta;        % x_dot depende de delta1, delta2
A_throttle(nx+1,nx+1) = -1/tau_throttle;     % delta1_dot
A_throttle(nx+2,nx+2) = -1/tau_throttle;     % delta2_dot

% Matriz B aumentada (mesmo número de entradas)
B_throttle = zeros(nx+2,nu);
B_throttle(1:nx,3:nu) = B_rest;      % efeito de i_t, delta_e, delta_a, delta_r
B_throttle(nx+1,1) = 1/tau_throttle; % u1 -> delta1_dot
B_throttle(nx+2,2) = 1/tau_throttle; % u2 -> delta2_dot

% Matriz C aumentada e D (mantendo saídas originais)
C_throttle = [C_base, zeros(ny,2)];
D_throttle = D_base;

% Autovalores sem dinâmicas de motores
lambda_sem = eig(A_base);

% Autovalores com dinâmicas de motores
lambda_com = eig(A_throttle);

fprintf('\nAutovalores da dinâmica linearizada SEM dinâmicas dos motores:\n\n');
damp(lambda_sem)

fprintf('\nAutovalores da dinâmica linearizada COM dinâmicas dos motores:\n\n');
damp(lambda_com)

% Procurar os dois polos adicionais devidos aos motores
% Esperamos aproximadamente -1/tau_throttle = -10
polos_reais = lambda_com(abs(imag(lambda_com)) < 1e-6);
polos_reais = sort(real(polos_reais));

fprintf('\nPolos reais (incluindo os dois associados às dinâmicas de throttle):\n');
for k = 1:length(polos_reais)
    fprintf('   % .4f\n', polos_reais(k));
end
fprintf('\nObserva-se a presença de dois polos próximos de -1/tau_throttle = -10 rad/s,\n');
fprintf('associados às dinâmicas de primeira ordem dos motores.\n');
fprintf('Os demais modos da aeronave praticamente não se alteram, como esperado.\n\n');

% Opcional: salvar o modelo aumentado
lin_output_L4_throttle = struct('A',A_throttle,'B',B_throttle, ...
                                'C',C_throttle,'D',D_throttle);
save('lin_L4_modelo_throttle.mat','lin_output_L4_throttle','X_eq','U_eq','Y_eq');
