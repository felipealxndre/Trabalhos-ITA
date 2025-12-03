clear
clc

% Gera modelo da Lista 3 (trim + linearização)
[A,B,C,D,X_eq,U_eq,Y_eq] = gera_modelo_L3();

fprintf('\n========== L4 - EXERCÍCIO 1: SAS LONGITUDINAL ==========\n');

% Índices dos estados: [1 V, 2 alpha, 3 q, 4 theta, 5 h, 6 x, 7 beta, 8 phi, 9 p, 10 r, 11 psi, 12 y]
i_alpha = 2;
i_q     = 3;

% Índice da entrada de profundor: u = [delta1; delta2; i_t; delta_e; delta_a; delta_r]
i_de = 4;

% Submodelo de período curto: x_sp = [alpha; q], u_sp = delta_e
A_sp = A([i_alpha i_q],[i_alpha i_q]);
B_sp = B([i_alpha i_q], i_de);

% Autovalores em malha aberta do período curto
lambda_sp_OL = eig(A_sp);
lam_OL = lambda_sp_OL(1);      % um dos dois (conjugados)

wn_sp_OL   = abs(lam_OL);
zeta_sp_OL = -real(lam_OL)/wn_sp_OL;

fprintf('\nPeríodo curto em malha aberta (modelo reduzido [alpha, q]):\n');
fprintf('   lambda_OL = % .4f %+.4fi\n', real(lambda_sp_OL(1)), imag(lambda_sp_OL(1)));
fprintf('                % .4f %+.4fi\n', real(lambda_sp_OL(2)), imag(lambda_sp_OL(2)));
fprintf('   wn_OL     = %.4f rad/s\n', wn_sp_OL);
fprintf('   zeta_OL   = %.4f\n\n', zeta_sp_OL);

% Especificações desejadas (Ex. 1(d) da AP3)
zeta_sp_des = 0.7;
wn_sp_des   = 3.0;

sigma_des = -zeta_sp_des*wn_sp_des;
wd_des    = wn_sp_des*sqrt(1 - zeta_sp_des^2);

lambda_sp_des = [sigma_des + 1i*wd_des, sigma_des - 1i*wd_des];

fprintf('Especificações desejadas para o período curto em malha fechada:\n');
fprintf('   zeta_des  = %.4f\n', zeta_sp_des);
fprintf('   wn_des    = %.4f rad/s\n', wn_sp_des);
fprintf('   lambda_des = % .4f %+.4fi\n', real(lambda_sp_des(1)), imag(lambda_sp_des(1)));
fprintf('                % .4f %+.4fi\n\n', real(lambda_sp_des(2)), imag(lambda_sp_des(2)));

% Projeto dos ganhos via alocação de polos (equivalente a escolher o ponto no RL)
K_sp = place(A_sp, B_sp, lambda_sp_des);
k_alpha = K_sp(1);
k_q     = K_sp(2);

fprintf('Ganhos do SAS longitudinal (lei de controle delta_e_SAS = -k_alpha*alpha - k_q*q):\n');
fprintf('   k_alpha = % .3f\n', k_alpha);
fprintf('   k_q     = % .3f\n\n', k_q);

% Matriz em malha fechada do modelo reduzido
A_sp_CL = A_sp - B_sp*K_sp;       % B_sp (2x1) * K_sp (1x2) -> 2x2

lambda_sp_CL = eig(A_sp_CL);
lam_CL = lambda_sp_CL(1);

wn_sp_CL   = abs(lam_CL);
zeta_sp_CL = -real(lam_CL)/wn_sp_CL;

fprintf('Período curto em malha fechada (modelo reduzido [alpha, q]):\n');
fprintf('   lambda_CL = % .4f %+.4fi\n', real(lambda_sp_CL(1)), imag(lambda_sp_CL(1)));
fprintf('                % .4f %+.4fi\n', real(lambda_sp_CL(2)), imag(lambda_sp_CL(2)));
fprintf('   wn_CL     = %.4f rad/s\n', wn_sp_CL);
fprintf('   zeta_CL   = %.4f\n\n', zeta_sp_CL);

% Aplicar os mesmos ganhos no modelo completo (12 estados)
K_full = zeros(1,size(A,1));
K_full(i_alpha) = k_alpha;
K_full(i_q)     = k_q;

A_full_CL = A - B(:,i_de)*K_full;     % realimentação apenas no canal do profundor

lambda_full_CL = eig(A_full_CL);

% Identificar o par de autovalores do período curto no modelo completo:
tol_im = 1e-6;
idx_cpx = find(abs(imag(lambda_full_CL)) > tol_im);
[~,idx_sort] = sort(abs(imag(lambda_full_CL(idx_cpx))),'descend');
idx_sp1 = idx_cpx(idx_sort(1));
idx_sp2 = idx_cpx(idx_sort(2));

lam1_full = lambda_full_CL(idx_sp1);
lam2_full = lambda_full_CL(idx_sp2);

wn_full_sp   = abs(lam1_full);
zeta_full_sp = -real(lam1_full)/wn_full_sp;

fprintf('Período curto em malha fechada (modelo completo 12 estados):\n');
fprintf('   lambda_CL_full = % .4f %+.4fi\n', real(lam1_full), imag(lam1_full));
fprintf('                     % .4f %+.4fi\n', real(lam2_full), imag(lam2_full));
fprintf('   wn_CL_full      = %.4f rad/s\n', wn_full_sp);
fprintf('   zeta_CL_full    = %.4f\n\n', zeta_full_sp);

fprintf('Conclusão: os ganhos posicionam o período curto próximo de zeta ≈ %.2f e wn ≈ %.2f rad/s,\n', ...
        zeta_sp_des, wn_sp_des);
fprintf('tanto no modelo reduzido quanto na dinâmica linearizada completa.\n\n');

% Salvar ganhos para uso no Ex.2
save('L4_Ex1_ganhos_SAS.mat','k_alpha','k_q');

