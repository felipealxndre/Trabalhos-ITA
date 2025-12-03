clear
close all
clc

load lin_output.mat

i_cond = 3;

% Full-order model (FOM):
A = lin_output(i_cond).A;
B = lin_output(i_cond).B;

% Reduced-order model (ROM):
A_sp = A(2:3,2:3);
B_sp = B(2:3,3);


%--------------------------------------------------------------------------
% 1A)

damp(A)
damp(A_sp)

% For the ROM:
% C_states = eye(size(A_sp,1));
% C_alpha = C_states(1,:);
C_alpha = [1 0];
ss_de_deg_2_alpha_deg = ss(A_sp,B_sp,C_alpha,0);
tf_de_deg_2_alpha_deg = tf(ss_de_deg_2_alpha_deg)
zpk_de_deg_2_alpha_deg = zpk(ss_de_deg_2_alpha_deg)

% % % For the FOM:
% C_states = eye(size(A,1));
% C_alpha = C_states(2,:);
% ss_de_deg_2_alpha_deg = ss(A,B(:,3),C_alpha,0);
% tf_de_deg_2_alpha_deg = tf(ss_de_deg_2_alpha_deg)
% zpk_de_deg_2_alpha_deg = zpk(ss_de_deg_2_alpha_deg)

% For the ROM:
% C_states = eye(size(A_sp,1));
% C_q = C_states(2,:);
C_q = [0 1];
ss_de_deg_2_q_deg_s = ss(A_sp,B_sp,C_q,0);
tf_de_deg_2_q_deg_s = tf(ss_de_deg_2_q_deg_s)
zpk_de_deg_2_q_deg_s = zpk(ss_de_deg_2_q_deg_s)

% % For the FOM:
% C_states = eye(size(A,1));
% C_q = C_states(3,:);
% ss_de_deg_2_q_deg_s = ss(A,B(:,3),C_q,0);
% tf_de_deg_2_q_deg_s = tf(ss_de_deg_2_q_deg_s)
% zpk_de_deg_2_q_deg_s = zpk(ss_de_deg_2_q_deg_s)


%--------------------------------------------------------------------------
% 1B) 

% k_vec = -logspace(-2,1,1000); %gera um vetor em que os valores são 10^-2 
% %até 1 com 1000 pontos
% % LEMBRAR QUE O GANHO É NEGATIVO
% 
% figure
% rlocus(ss_de_deg_2_alpha_deg,k_vec)
% 
% 
% %--------------------------------------------------------------------------
% % 1C)
% 
% k_vec = -logspace(-2,1,1000);
% 
% figure
% rlocus(ss_de_deg_2_q_deg_s,k_vec)


%--------------------------------------------------------------------------
% 1D)

k_vec = -logspace(-2,1,1000);

figure
rlocus(ss_de_deg_2_alpha_deg,k_vec)

k_alpha = -2.33; %escolheu o ganho que deu 10% a menos que o requisito de 
%wn, portanto que deu wn = 2.7 (ganho sempre negativo)

% ss_de_deg_2_alpha_deg = ss(A_sp,B_sp,C_alpha,0);
A_cl_alpha = A_sp - B_sp*k_alpha*C_alpha;

% Root locus for feedback of q to delta_e after loop closure for alpha:
ss_de_deg_2_q_deg_s = ss(A_cl_alpha,B_sp,C_q,0);

k_vec = -logspace(-2,1,1000);

figure
rlocus(ss_de_deg_2_q_deg_s,k_vec)

k_q = -1.35;

A_sp_cl = A_cl_alpha - B_sp*k_q*C_q;

OHara_plot

[wn_sp_ol,zeta_sp_ol] = damp(A_sp);
plot(zeta_sp_ol(2),wn_sp_ol(2),'wx','MarkerSize',8,'LineWidth',2)

[wn_cl_alpha,zeta_cl_alpha] = damp(A_cl_alpha);
plot(zeta_cl_alpha(2),wn_cl_alpha(2),'rx','MarkerSize',9,'LineWidth',2)

[wn_sp_cl,zeta_sp_cl] = damp(A_sp_cl);
plot(zeta_sp_cl(2),wn_sp_cl(2),'bx','MarkerSize',10,'LineWidth',2)


%--------------------------------------------------------------------------
% 1E)

K = [k_alpha k_q];

save K_1D.mat K

C_states = eye(size(A,1));
C_alpha = C_states(2,:);
C_q = C_states(3,:);

C = [
    C_alpha
    C_q
    ];

A_cl = A - B(:, 3)*K*C;

damp(A)
damp(A_cl)
damp(A_sp_cl)

sel_states = [1:5 7:11];

[eigvec,eigval] = eig(A(sel_states,sel_states));
eigval = diag(eigval);
[eigval,i_sort] = sort(eigval);
eigvec = eigvec(:,i_sort);
damp(eigval)

[eigvec_cl,eigval_cl] = eig(A_cl(sel_states,sel_states));
eigval_cl = diag(eigval_cl);
[eigval_cl,i_sort] = sort(eigval_cl);
eigvec_cl = eigvec_cl(:,i_sort);
damp(eigval_cl)

damp(A_sp_cl)


%--------------------------------------------------------------------------
% 1F)


%vetor de saída com 5 componentes. delta_alpha, delta_q,
%deltadelta_e_piloto, deltadelta_e_SAS, deltadelta_e

% y = [alpha_deg q_deg_s (delta_e_deg)_pilot (delta_e_deg)_SAS delta_e_deg].';



C_SAS = C;

% Feedback of alpha and q:
C = [
    C_alpha
    C_q
    zeros(1,size(A,1))
    -K*C_SAS
    -K*C_SAS
    ];
D = [
    zeros(1,size(B,2))
    zeros(1,size(B,2))
    0 0 1 0 0
    zeros(1,size(B,2))
    0 0 1 0 0
    ];

ss_sim_cl = ss(A_cl,B,C,D);

% Open loop:
K_sim = zeros(1,2); %no open loop basta zerar os ganhos!
C = [
    C_alpha
    C_q
    zeros(1,size(A,1))
    -K_sim*C_SAS
    -K_sim*C_SAS
    ];
D = [
    zeros(1,size(B,2))
    zeros(1,size(B,2))
    0 0 1 0 0
    zeros(1,size(B,2))
    0 0 1 0 0
    ];

ss_sim_ol = ss(A,B,C,D);

% Feedback of alpha only:
K_sim = [k_alpha 0]; %ganho de alpha apenas então zera ganho de k_q
C = [
    C_alpha
    C_q
    zeros(1,size(A,1))
    -K_sim*C_SAS
    -K_sim*C_SAS
    ];
D = [
    zeros(1,size(B,2))
    zeros(1,size(B,2))
    0 0 1 0 0
    zeros(1,size(B,2))
    0 0 1 0 0
    ];

A_cl_alpha = A - B(:,3)*K_sim*C_SAS;

ss_sim_cl_alpha = ss(A_cl_alpha,B,C,D);

% Feedback of q only:
K_sim = [0 k_q]; %ganho de k_alpha zerado
C = [
    C_alpha
    C_q
    zeros(1,size(A,1))
    -K_sim*C_SAS
    -K_sim*C_SAS
    ];
D = [
    zeros(1,size(B,2))
    zeros(1,size(B,2))
    0 0 1 0 0
    zeros(1,size(B,2))
    0 0 1 0 0
    ];

A_cl_q = A - B(:,3)*K_sim*C_SAS;

ss_sim_cl_q = ss(A_cl_q,B,C,D);

% Simulation of initial condition:
x0 = zeros(size(A,1),1);
x0(2) = 1; %perturbação de 1 grau em alpha

t_f = 5;
dt = 0.010;

[y_cl,t_cl,x_cl] = initial(ss_sim_cl,x0,0:dt:t_f);
[y_ol,t_ol,x_ol] = initial(ss_sim_ol,x0,0:dt:t_f);
[y_cl_alpha,t_cl_alpha,x_cl_alpha] = initial(ss_sim_cl_alpha,x0,0:dt:t_f);
[y_cl_q,t_cl_q,x_cl_q] = initial(ss_sim_cl_q,x0,0:dt:t_f);

figure
subplot(311)
plot(t_ol,y_ol(:,1))
hold all
plot(t_cl_alpha,y_cl_alpha(:,1))
plot(t_cl_q,y_cl_q(:,1))
plot(t_cl,y_cl(:,1))
grid on
xlabel('t [s]')
ylabel('\Delta\alpha [deg]')

subplot(312)
plot(t_ol,y_ol(:,2))
hold all
plot(t_cl_alpha,y_cl_alpha(:,2))
plot(t_cl_q,y_cl_q(:,2))
plot(t_cl,y_cl(:,2))
grid on
xlabel('t [s]')
ylabel('q [deg/s]')

subplot(313)
plot(t_ol,y_ol(:,5))
hold all
plot(t_cl_alpha,y_cl_alpha(:,5))
plot(t_cl_q,y_cl_q(:,5))
plot(t_cl,y_cl(:,5))
grid on
xlabel('t [s]')
ylabel('\Delta\delta_e [deg]')

legend('OL','k_\alpha only','k_q only','CL','Location','NorthEast','Orientation','Horizontal')


%--------------------------------------------------------------------------
% 1G)

% Simulation of initial condition:
x0 = zeros(size(A,1),1);
x0(3) = 1;

t_f = 5;
dt = 0.010;

[y_cl,t_cl,x_cl] = initial(ss_sim_cl,x0,0:dt:t_f);
[y_ol,t_ol,x_ol] = initial(ss_sim_ol,x0,0:dt:t_f);
[y_cl_alpha,t_cl_alpha,x_cl_alpha] = initial(ss_sim_cl_alpha,x0,0:dt:t_f);
[y_cl_q,t_cl_q,x_cl_q] = initial(ss_sim_cl_q,x0,0:dt:t_f);

figure
subplot(311)
plot(t_ol,y_ol(:,1))
hold all
plot(t_cl_alpha,y_cl_alpha(:,1))
plot(t_cl_q,y_cl_q(:,1))
plot(t_cl,y_cl(:,1))
grid on
xlabel('t [s]')
ylabel('\Delta\alpha [deg]')

subplot(312)
plot(t_ol,y_ol(:,2))
hold all
plot(t_cl_alpha,y_cl_alpha(:,2))
plot(t_cl_q,y_cl_q(:,2))
plot(t_cl,y_cl(:,2))
grid on
xlabel('t [s]')
ylabel('q [deg/s]')

subplot(313)
plot(t_ol,y_ol(:,5))
hold all
plot(t_cl_alpha,y_cl_alpha(:,5))
plot(t_cl_q,y_cl_q(:,5))
plot(t_cl,y_cl(:,5))
grid on
xlabel('t [s]')
ylabel('\Delta\delta_e [deg]')

legend('OL','k_\alpha only','k_q only','CL','Location','NorthEast','Orientation','Horizontal')


%--------------------------------------------------------------------------
% 1H)

dt = 0.010;
t_f = 20;
t = 0:dt:t_f;

%comando doublet
delta_e_deg_vec = zeros(length(t),1);
delta_e_deg_vec(t>=1) = 3;
delta_e_deg_vec(t>=3) = -3;
delta_e_deg_vec(t>=5) = 0;

x0 = zeros(size(A,1),1);

%simulação de resposta da dinâmica linear e não é mais initial 
[y_cl,t_cl,x_cl] = lsim(ss_sim_cl(:,3),delta_e_deg_vec,t,x0);
[y_ol,t_ol,x_ol] = lsim(ss_sim_ol(:,3),delta_e_deg_vec,t,x0);
[y_cl_alpha,t_cl_alpha,x_cl_alpha] = lsim(ss_sim_cl_alpha(:,3),delta_e_deg_vec,t,x0);
[y_cl_q,t_cl_q,x_cl_q] = lsim(ss_sim_cl_q(:,3),delta_e_deg_vec,t,x0);

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
plot(t_ol,y_ol(:,5))
hold all
plot(t_cl_alpha,y_cl_alpha(:,5))
plot(t_cl_q,y_cl_q(:,5))
plot(t_cl,y_cl(:,5))
grid on
xlabel('t [s]')
ylabel('\Delta\delta_e [deg]')

subplot(514)
plot(t_ol,y_ol(:,3))
hold all
plot(t_cl_alpha,y_cl_alpha(:,3))
plot(t_cl_q,y_cl_q(:,3))
plot(t_cl,y_cl(:,3))
grid on
xlabel('t [s]')
ylabel('(\Delta\delta_e)_{pilot} [deg]')

subplot(515)
plot(t_ol,y_ol(:,4))
hold all
plot(t_cl_alpha,y_cl_alpha(:,4))
plot(t_cl_q,y_cl_q(:,4))
plot(t_cl,y_cl(:,4))
grid on
xlabel('t [s]')
ylabel('(\Delta\delta_e)_{SAS} [deg]')

legend('OL','k_\alpha only','k_q only','CL','Location','NorthEast','Orientation','Horizontal')

