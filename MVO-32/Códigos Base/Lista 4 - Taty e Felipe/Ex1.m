clear
close all
clc

[A,B,C,D,X_eq,U_eq,Y_eq] = gera_modelo_L3();

A_sp = A(2:3,2:3);
B_sp = B(2:3,4);

C_alpha = [1 0];
C_q     = [0 1];

fprintf('\n=== MALHA ABERTA (MODELO REDUZIDO) ===\n');
lambda_ol = eig(A_sp);
fprintf('lambda_1 = % .4f %+.4fi\n', real(lambda_ol(1)), imag(lambda_ol(1)));
fprintf('lambda_2 = % .4f %+.4fi\n', real(lambda_ol(2)), imag(lambda_ol(2)));

[wn_sp_ol,zeta_sp_ol] = damp(A_sp);
fprintf('wn_OL   = %.4f  %.4f rad/s\n', wn_sp_ol(1), wn_sp_ol(2));
fprintf('zeta_OL = %.4f  %.4f\n', zeta_sp_ol(1), zeta_sp_ol(2));


k_vec = -logspace(-2,1,1000);
ss_de_deg_2_alpha_deg = ss(A_sp,B_sp,C_alpha,0);

figure
rlocus(ss_de_deg_2_alpha_deg,k_vec)

k_alpha = -0.336;

A_cl_alpha = A_sp - B_sp*k_alpha*C_alpha;

fprintf('\n=== APÓS FECHAR LOOP EM ALPHA ===\n');
lambda_alpha = eig(A_cl_alpha);
fprintf('lambda_1 = % .4f %+.4fi\n', real(lambda_alpha(1)), imag(lambda_alpha(1)));
fprintf('lambda_2 = % .4f %+.4fi\n', real(lambda_alpha(2)), imag(lambda_alpha(2)));

[wn_cl_alpha,zeta_cl_alpha] = damp(A_cl_alpha);
fprintf('wn_alpha_CL   = %.4f  %.4f rad/s\n', wn_cl_alpha(1), wn_cl_alpha(2));
fprintf('zeta_alpha_CL = %.4f  %.4f\n', zeta_cl_alpha(1), zeta_cl_alpha(2));


ss_de_deg_2_q_deg_s = ss(A_cl_alpha,B_sp,C_q,0);
k_vec = -logspace(-2,1,1000);

figure
rlocus(ss_de_deg_2_q_deg_s,k_vec)

k_q = -0.540;

A_sp_cl = A_cl_alpha - B_sp*k_q*C_q;

fprintf('\n=== MALHA FECHADA FINAL (ALPHA + Q) ===\n');
lambda_cl = eig(A_sp_cl);
fprintf('lambda_1 = % .4f %+.4fi\n', real(lambda_cl(1)), imag(lambda_cl(1)));
fprintf('lambda_2 = % .4f %+.4fi\n', real(lambda_cl(2)), imag(lambda_cl(2)));

[wn_sp_cl,zeta_sp_cl] = damp(A_sp_cl);
fprintf('wn_CL   = %.4f  %.4f rad/s\n', wn_sp_cl(1), wn_sp_cl(2));
fprintf('zeta_CL = %.4f  %.4f\n', zeta_sp_cl(1), zeta_sp_cl(2));

K_full = zeros(1,12);
K_full(2) = -0.336;   % k_alpha
K_full(3) = -0.540;   % k_q

A_full_CL = A - B(:,4)*K_full;

lambda_full = eig(A_full_CL);

fprintf('\n=== AUTOVALORES DO PERÍODO CURTO (MODELO COMPLETO) ===\n');
for i = 1:12
    fprintf('lambda_%d = % .4f %+.4fi\n', i, real(lambda_full(i)), imag(lambda_full(i)));
end

OHara_plot
plot(zeta_sp_ol(2),wn_sp_ol(2),'wx','MarkerSize',8,'LineWidth',2)
plot(zeta_cl_alpha(2),wn_cl_alpha(2),'rx','MarkerSize',9,'LineWidth',2)
plot(zeta_sp_cl(2),wn_sp_cl(2),'bx','MarkerSize',10,'LineWidth',2)

fprintf('\nGanhos finais:\n');
fprintf('k_alpha = % .4f\n', k_alpha);
fprintf('k_q     = % .4f\n\n', k_q);

