clear
close all
clc

[A,B,C,D,X_eq,U_eq,Y_eq] = gera_modelo_L3();
load('L4_Ex1_ganhos_SAS.mat','k_alpha','k_q')

i_alpha = 2;
i_q     = 3;
i_de    = 4;

n = size(A,1);
rad2deg = 180/pi;

B_de_deg = B(:,i_de)*pi/180;

e_alpha = zeros(1,n); e_alpha(i_alpha) = 1;
e_q     = zeros(1,n); e_q(i_q)       = 1;

K_OL    = zeros(1,n);
K_alpha = zeros(1,n); K_alpha(i_alpha) = k_alpha;
K_q     = zeros(1,n); K_q(i_q)         = k_q;
K_both  = zeros(1,n); K_both(i_alpha)  = k_alpha; K_both(i_q) = k_q;

A_OL    = A;
A_alpha = A - B_de_deg*K_alpha;
A_q     = A - B_de_deg*K_q;
A_both  = A - B_de_deg*K_both;

C_OL = [rad2deg*e_alpha;
        rad2deg*e_q;
        -rad2deg*K_OL;
        zeros(1,n);
        -rad2deg*K_OL];
D_OL = [0;0;1;1;0];

C_alpha_sys = [rad2deg*e_alpha;
               rad2deg*e_q;
               -rad2deg*K_alpha;
               zeros(1,n);
               -rad2deg*K_alpha];
D_alpha_sys = [0;0;1;1;0];

C_q_sys = [rad2deg*e_alpha;
           rad2deg*e_q;
           -rad2deg*K_q;
           zeros(1,n);
           -rad2deg*K_q];
D_q_sys = [0;0;1;1;0];

C_both_sys = [rad2deg*e_alpha;
              rad2deg*e_q;
              -rad2deg*K_both;
              zeros(1,n);
              -rad2deg*K_both];
D_both_sys = [0;0;1;1;0];

ss_sim_ol       = ss(A_OL,   B_de_deg,C_OL,       D_OL);
ss_sim_cl_alpha = ss(A_alpha,B_de_deg,C_alpha_sys,D_alpha_sys);
ss_sim_cl_q     = ss(A_q,    B_de_deg,C_q_sys,    D_q_sys);
ss_sim_cl       = ss(A_both, B_de_deg,C_both_sys, D_both_sys);

dt = 0.010;
t_f = 20;
t = 0:dt:t_f;

delta_e_deg_vec = zeros(length(t),1);
delta_e_deg_vec(t>=1) = 3;
delta_e_deg_vec(t>=3) = -3;
delta_e_deg_vec(t>=5) = 0;

x0 = zeros(n,1);

[y_cl,t_cl,x_cl]             = lsim(ss_sim_cl,      delta_e_deg_vec,t,x0);
[y_ol,t_ol,x_ol]             = lsim(ss_sim_ol,      delta_e_deg_vec,t,x0);
[y_cl_alpha,t_cl_alpha,x_cl_alpha] = lsim(ss_sim_cl_alpha,delta_e_deg_vec,t,x0);
[y_cl_q,t_cl_q,x_cl_q]       = lsim(ss_sim_cl_q,    delta_e_deg_vec,t,x0);

fprintf('\nMáximos de |Δα| e |q| (em graus) para doublet de 3° em δ_e:\n\n');

fprintf('Malha aberta:\n');
fprintf('   max |Δα| = %6.3f deg\n',   max(abs(y_ol(:,1))));
fprintf('   max |q|  = %6.3f deg/s\n\n', max(abs(y_ol(:,2))));

fprintf('SAS apenas em α:\n');
fprintf('   max |Δα| = %6.3f deg\n',   max(abs(y_cl_alpha(:,1))));
fprintf('   max |q|  = %6.3f deg/s\n\n', max(abs(y_cl_alpha(:,2))));

fprintf('SAS apenas em q:\n');
fprintf('   max |Δα| = %6.3f deg\n',   max(abs(y_cl_q(:,1))));
fprintf('   max |q|  = %6.3f deg/s\n\n', max(abs(y_cl_q(:,2))));

fprintf('SAS em α e q:\n');
fprintf('   max |Δα| = %6.3f deg\n',   max(abs(y_cl(:,1))));
fprintf('   max |q|  = %6.3f deg/s\n\n', max(abs(y_cl(:,2))));

figure
subplot(511)
plot(t_ol,y_ol(:,1))
hold all
plot(t_cl_alpha,y_cl_alpha(:,1))
plot(t_cl_q,y_cl_q(:,1))
plot(t_cl,y_cl(:,1))
grid on
xlabel('t [s]')
ylabel('\Delta\alpha [deg]')

subplot(512)
plot(t_ol,y_ol(:,2))
hold all
plot(t_cl_alpha,y_cl_alpha(:,2))
plot(t_cl_q,y_cl_q(:,2))
plot(t_cl,y_cl(:,2))
grid on
xlabel('t [s]')
ylabel('q [deg/s]')

subplot(513)
plot(t_ol,y_ol(:,3))
hold all
plot(t_cl_alpha,y_cl_alpha(:,3))
plot(t_cl_q,y_cl_q(:,3))
plot(t_cl,y_cl(:,3))
grid on
xlabel('t [s]')
ylabel('\Delta\delta_e [deg]')

subplot(514)
plot(t_ol,y_ol(:,4))
hold all
plot(t_cl_alpha,y_cl_alpha(:,4))
plot(t_cl_q,y_cl_q(:,4))
plot(t_cl,y_cl(:,4))
grid on
xlabel('t [s]')
ylabel('(\Delta\delta_e)_{pilot} [deg]')

subplot(515)
plot(t_ol,y_ol(:,5))
hold all
plot(t_cl_alpha,y_cl_alpha(:,5))
plot(t_cl_q,y_cl_q(:,5))
plot(t_cl,y_cl(:,5))
grid on
xlabel('t [s]')
ylabel('(\Delta\delta_e)_{SAS} [deg]')

legend('OL','k_\alpha only','k_q only','CL','Location','NorthEast','Orientation','Horizontal')
