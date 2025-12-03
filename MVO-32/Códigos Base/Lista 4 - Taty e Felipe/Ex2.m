clear
clc

% Gera modelo da Lista 3
[A,B,C,D,X_eq,U_eq,Y_eq] = gera_modelo_L3();

% Carrega ganhos do SAS projetados no Ex.1
load('L4_Ex1_ganhos_SAS.mat','k_alpha','k_q');

fprintf('\n========== L4 - EXERCÍCIO 2: REPRODUÇÃO DO EX. 1(h) DA AP3 ==========\n');

% Índices de estados e entrada
i_alpha = 2;
i_q     = 3;
i_de    = 4;      % coluna de delta_e em B

B_de = B(:,i_de);

% Tempo e entrada do piloto (doublet de 3 deg)
dt = 0.01;
t_final = 20;
t = (0:dt:t_final)';

delta_deg = 3;             % 3 graus
delta_rad = delta_deg*pi/180;

u_pilot = zeros(length(t),1);
u_pilot(t>=1 & t<3) =  delta_rad;
u_pilot(t>=3 & t<5) = -delta_rad;

% Caso 0: malha aberta (sem SAS)
A_OL = A;
B_OL = B_de;

sys_OL = ss(A_OL,B_OL,eye(size(A_OL)),zeros(size(A_OL,1),1));
[~,~,x_OL] = lsim(sys_OL,u_pilot,t);

alpha_OL = x_OL(:,i_alpha);    % rad
q_OL     = x_OL(:,i_q);        % rad/s

% Caso 1: SAS apenas em alpha (k_q = 0)
K_full_alpha = zeros(1,size(A,1));
K_full_alpha(i_alpha) = k_alpha;

A_CL_alpha = A - B_de*K_full_alpha;
B_CL_alpha = B_de;    % entrada ainda é o delta_e_pilot

sys_CL_alpha = ss(A_CL_alpha,B_CL_alpha,eye(size(A,1)),zeros(size(A,1),1));
[~,~,x_CL_alpha] = lsim(sys_CL_alpha,u_pilot,t);

alpha_CL_alpha = x_CL_alpha(:,i_alpha);
q_CL_alpha     = x_CL_alpha(:,i_q);
de_SAS_alpha   = -(K_full_alpha * x_CL_alpha.').';   % saída SAS
de_total_alpha = u_pilot + de_SAS_alpha;

% Caso 2: SAS apenas em q (k_alpha = 0)
K_full_q = zeros(1,size(A,1));
K_full_q(i_q) = k_q;

A_CL_q = A - B_de*K_full_q;
B_CL_q = B_de;

sys_CL_q = ss(A_CL_q,B_CL_q,eye(size(A,1)),zeros(size(A,1),1));
[~,~,x_CL_q] = lsim(sys_CL_q,u_pilot,t);

alpha_CL_q = x_CL_q(:,i_alpha);
q_CL_q     = x_CL_q(:,i_q);
de_SAS_q   = -(K_full_q * x_CL_q.').';
de_total_q = u_pilot + de_SAS_q;

% Caso 3: SAS em alpha e q
K_full_both = zeros(1,size(A,1));
K_full_both(i_alpha) = k_alpha;
K_full_both(i_q)     = k_q;

A_CL_both = A - B_de*K_full_both;
B_CL_both = B_de;

sys_CL_both = ss(A_CL_both,B_CL_both,eye(size(A,1)),zeros(size(A,1),1));
[~,~,x_CL_both] = lsim(sys_CL_both,u_pilot,t);

alpha_CL_both = x_CL_both(:,i_alpha);
q_CL_both     = x_CL_both(:,i_q);
de_SAS_both   = -(K_full_both * x_CL_both.').';
de_total_both = u_pilot + de_SAS_both;

% Conversão para graus
rad2deg = 180/pi;

alpha_OL_deg        = alpha_OL*rad2deg;
q_OL_deg            = q_OL*rad2deg;
alpha_CL_alpha_deg  = alpha_CL_alpha*rad2deg;
q_CL_alpha_deg      = q_CL_alpha*rad2deg;
alpha_CL_q_deg      = alpha_CL_q*rad2deg;
q_CL_q_deg          = q_CL_q*rad2deg;
alpha_CL_both_deg   = alpha_CL_both*rad2deg;
q_CL_both_deg       = q_CL_both*rad2deg;
de_pilot_deg        = u_pilot*rad2deg;
de_SAS_alpha_deg    = de_SAS_alpha*rad2deg;
de_SAS_q_deg        = de_SAS_q*rad2deg;
de_SAS_both_deg     = de_SAS_both*rad2deg;
de_total_alpha_deg  = de_total_alpha*rad2deg;
de_total_q_deg      = de_total_q*rad2deg;
de_total_both_deg   = de_total_both*rad2deg;

% Impressão de alguns valores numéricos (máximos)

fprintf('\nMáximos de |alpha| e |q| (em graus) para o doublet de 3° no profundor:\n\n');

fprintf('Malha aberta:\n');
fprintf('   max |alpha| = %6.3f deg\n', max(abs(alpha_OL_deg)));
fprintf('   max |q|     = %6.3f deg/s\n\n', max(abs(q_OL_deg)));

fprintf('SAS apenas em alpha:\n');
fprintf('   max |alpha| = %6.3f deg\n', max(abs(alpha_CL_alpha_deg)));
fprintf('   max |q|     = %6.3f deg/s\n\n', max(abs(q_CL_alpha_deg)));

fprintf('SAS apenas em q:\n');
fprintf('   max |alpha| = %6.3f deg\n', max(abs(alpha_CL_q_deg)));
fprintf('   max |q|     = %6.3f deg/s\n\n', max(abs(q_CL_q_deg)));

fprintf('SAS em alpha e q:\n');
fprintf('   max |alpha| = %6.3f deg\n', max(abs(alpha_CL_both_deg)));
fprintf('   max |q|     = %6.3f deg/s\n\n', max(abs(q_CL_both_deg)));

fprintf('Resultados numéricos prontos. Gerando gráficos para o relatório...\n');

%% GRÁFICOS – REPRODUÇÃO DO EX. 1(h) DA AP3

% 1) alpha(t) em graus – quatro configurações
figure;
plot(t, alpha_OL_deg,        'k-',  ...  % malha aberta
     t, alpha_CL_alpha_deg,  'b--', ...
     t, alpha_CL_q_deg,      'r-.', ...
     t, alpha_CL_both_deg,   'g:');
grid on;
xlabel('t [s]');
ylabel('\alpha [deg]');
title('\alpha(t) para doublet de 3^\circ em \delta_e');
legend('Malha aberta', 'SAS em \alpha', 'SAS em q', 'SAS em \alpha e q', ...
       'Location','Best');

% 2) q(t) em graus/s – quatro configurações
figure;
plot(t, q_OL_deg,        'k-',  ...
     t, q_CL_alpha_deg,  'b--', ...
     t, q_CL_q_deg,      'r-.', ...
     t, q_CL_both_deg,   'g:');
grid on;
xlabel('t [s]');
ylabel('q [deg/s]');
title('q(t) para doublet de 3^\circ em \delta_e');
legend('Malha aberta', 'SAS em \alpha', 'SAS em q', 'SAS em \alpha e q', ...
       'Location','Best');

% 3) Comandos de profundor – piloto, SAS e total (caso SAS em \alpha e q)
figure;
plot(t, de_pilot_deg,      'k-',  ...
     t, de_SAS_both_deg,   'r--', ...
     t, de_total_both_deg, 'b-.');
grid on;
xlabel('t [s]');
ylabel('\delta_e [deg]');
title('Comandos de profundor: piloto, SAS e total');
legend('\delta_{e,pilot}', '(\delta_e)_{SAS}', '\delta_{e,total}', ...
       'Location','Best');