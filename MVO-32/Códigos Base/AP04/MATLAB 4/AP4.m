clear
close all
clc

load lin_output.mat

i_cond = 3;

A_ld = lin_output(i_cond).A(7:10,7:10);
B_ld = lin_output(i_cond).B(7:10,4:5);

n = size(A_ld,1);
m = size(B_ld,2);

%--------------------------------------------------------------------------
% 1A)

tau_a = 0.050;
tau_r = 0.050;
tau_w = 1;

C_states = eye(n);
C_r = C_states(4,:);

A = [
    
    ];
B = [
    
    ];

[eigvec_ld,eigval_ld] = eig(A_ld);
eigval_ld = diag(eigval_ld);

[eigvec,eigval] = eig(A);
eigval = diag(eigval);

damp(eigval_ld)
damp(eigval)

%--------------------------------------------------------------------------
% 1B)

Q = diag([1 2 2 1 0 0 0]);
R = 2*eye(2);

eig(Q)
eig(R)

% z = [beta phi p r].';

Qbar = zeros(4,4);

Qbar(1,1) = 1/3^2;
Qbar(2,2) = 1/5^2;
Qbar(3,3) = 1/5^2;
Qbar(4,4) = 1/3^2;

C_states = eye(size(A,1));
C_beta = C_states(1,:);
C_phi = C_states(2,:);
C_p = C_states(3,:);
C_r = C_states(4,:);
C_x_w = C_states(7,:);

% z = H*x;

H = [
    C_beta
    C_phi
    C_p
    C_r
    ];

Q = H.'*Qbar*H

R(1,1) = 1/2^2;
R(2,2) = 1/5^2;

R

eig(Q)
eig(R)

%--------------------------------------------------------------------------
% 1C)

fK = @(k)([]);

k0 = 

K0 = fK(k0)

C_r_w = C_x_w + C_r;

C = [
    C_p
    C_r_w
    ];

damp(A - B*K0*C)

%--------------------------------------------------------------------------
% 1D)

fPI = @(k)(PI(fK(k),A,B,C,Q,R));

options = optimset('Display','iter','TolX',1e-6,'TolFun',1e-6);

k = fminsearch(fPI,k0,options);

K = fK(k)

A_c = A - B*K*C;

damp(A_c)

figure
plot(real(eig(A)),imag(eig(A)),'bx','MarkerSize',8,'LineWidth',2)
hold on
grid on
plot(real(eig(A_c)),imag(eig(A_c)),'rx','MarkerSize',7,'LineWidth',2)
xlabel('Real part')
ylabel('Imaginary part')
legend('OL','CL')

[eigvec_c,eigval_c] = eig(A_c);
eigval_c = diag(eigval_c);

damp(eigval_c)

Jf = fPI(k)

fPI_Q = @(k)(PI(fK(k),A,B,C,Q,0*R));
Jf_Q = fPI_Q(k)

fPI_R = @(k)(PI(fK(k),A,B,C,0*Q,R));
Jf_R = fPI_R(k)

Jf_Q+Jf_R

K_1D = K;
C_1D = C;

Jf_1D = Jf;

%--------------------------------------------------------------------------
% 1E)

fK_1E = @(k)([]);

C_1E = [
    C_phi
    C_p
    C_beta
    C_r_w
    ];

fPI_1E = @(k)(PI(fK_1E(k),A,B,C_1E,Q,R));

k0_1E = ;

k_1E = fminsearch(fPI_1E,k0_1E,options);

K_1E = fK_1E(k_1E)

A_c_1E = A - B*K_1E*C_1E;

damp(A_c_1E)

figure
plot(real(eig(A)),imag(eig(A)),'bx','MarkerSize',8,'LineWidth',2)
hold on
grid on
plot(real(eig(A_c)),imag(eig(A_c)),'rx','MarkerSize',7,'LineWidth',2)
plot(real(eig(A_c_1E)),imag(eig(A_c_1E)),'kx','MarkerSize',6,'LineWidth',2)
xlabel('Real part')
ylabel('Imaginary part')
legend('OL','CL 1D','CL 1E')

[eigvec_c_1E,eigval_c_1E] = eig(A_c_1E);
eigval_c_1E = diag(eigval_c_1E);

damp(eigval_c_1E)

Jf_1E = fPI_1E(k_1E);

%--------------------------------------------------------------------------
% 1F)

Sim_Plot_1F

%--------------------------------------------------------------------------
% 1G)

Jf_1D
Jf_1E
