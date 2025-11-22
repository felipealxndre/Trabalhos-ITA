clear
close all
clc

global g
global aircraft
create_aircraft;

g = 9.80665;

trim_par = struct('V',228.8138886191133,'h',11000);

options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);

%% Ex.2

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

%% Ex.3

nX_Ex3 = length(X_eq_Ex2);
nU_Ex3 = length(U_eq_Ex2);
nY_Ex3 = length(Y_eq_Ex2);

lin_output_Ex3 = struct('A',zeros(nX_Ex3,nX_Ex3),'B',zeros(nX_Ex3,nU_Ex3),'C',zeros(nY_Ex3,nX_Ex3),'D',zeros(nY_Ex3,nU_Ex3));

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

lin_output_Ex3.A = A_Ex3; lin_output_Ex3.B = B_Ex3; lin_output_Ex3.C = C_Ex3; lin_output_Ex3.D = D_Ex3;

sel_states_Ex3 = 1:12;
[eigvec_Ex3,eigval_Ex3] = eig(lin_output_Ex3.A(sel_states_Ex3,sel_states_Ex3));
eigval_Ex3 = diag(eigval_Ex3);
[eigval_Ex3,i_sort_Ex3] = sort(eigval_Ex3);
eigvec_Ex3 = eigvec_Ex3(:,i_sort_Ex3);

fprintf('\n----- EX.3 -----\n\n');
damp(eigval_Ex3)

%% Ex.4

<<<<<<< HEAD:MVO-32/Códigos Base/Lista 3 - Taty e Felipe/MAIN.m
x0_Ex4 = [ X_eq_Ex2(2); U_eq_Ex2(1); U_eq_Ex2(3); 0.0; 0.5; -1.0];
=======
x0_Ex4 = [
    X_eq_Ex2(1);    % V
    X_eq_Ex2(2);    % alpha
    0;              % q
    X_eq_Ex2(4);    % theta
    0;              % phi
    0;              % p
    0;              % r
    0;              % psi
    U_eq_Ex2(1);    % delta_1
    U_eq_Ex2(3);    % i_t
    0;              % delta_a
    0               % delta_r
];

% Aumentar 'MaxFunEvals' e 'MaxIter' para dar fôlego ao otimizador 12x12
options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10, 'MaxFunEvals', 10000, 'MaxIter', 1000);

>>>>>>> 6ce0615 (feat: Relatório de Ensaio em Voo):MVO-32/Códigos Base/MAIN.m
[x_sol_Ex4, ~] = fsolve(@trim_function_Ex4, x0_Ex4, options, trim_par);
[~, X_eq_Ex4, U_eq_Ex4, Y_eq_Ex4] = trim_function_Ex4(x_sol_Ex4, trim_par);
fprintf('\n----- EX.4 -----\n\n');
fprintf('   %-12s = %8.2f %s\n','V',X_eq_Ex4(1),'m/s');
fprintf('   %-12s = %8.3f %s\n','alpha',X_eq_Ex4(2),'deg');
fprintf('   %-12s = %8.3f %s\n','q',X_eq_Ex4(3),'deg/s');
fprintf('   %-12s = %8.3f %s\n','theta',X_eq_Ex4(4),'deg');
fprintf('   %-12s = %8.1f %s\n','h',X_eq_Ex4(5),'m');
fprintf('   %-12s = %8.3f %s\n','beta',X_eq_Ex4(7),'deg');
fprintf('   %-12s = %8.3f %s\n','phi',X_eq_Ex4(8),'deg');
fprintf('   %-12s = %8.3f %s\n','p',X_eq_Ex4(9),'deg/s');
fprintf('   %-12s = %8.3f %s\n','r',X_eq_Ex4(10),'deg/s');
fprintf('\n');
fprintf('   %-12s = %8.2f %s\n','delta1',U_eq_Ex4(1),'');
fprintf('   %-12s = %8.2f %s\n','delta2',U_eq_Ex4(2),'');
fprintf('   %-12s = %8.3f %s\n','i_t',U_eq_Ex4(3),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_e',U_eq_Ex4(4),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_a',U_eq_Ex4(5),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_r',U_eq_Ex4(6),'deg');
fprintf('\n');
fprintf('   %-12s = %8.3f %s\n','gamma',Y_eq_Ex4(1),'deg');
fprintf('   %-12s = %8.1f %s\n','T1',Y_eq_Ex4(2),'N');
fprintf('   %-12s = %8.1f %s\n','T2',Y_eq_Ex4(3),'N');
fprintf('   %-12s = %8.3f %s\n','Mach',Y_eq_Ex4(4),'');
fprintf('   %-12s = %8.4f %s\n','C_D',Y_eq_Ex4(5),'');
fprintf('   %-12s = %8.4f %s\n','C_L',Y_eq_Ex4(6),'');
fprintf('   %-12s = %8.4f %s\n','C_m',Y_eq_Ex4(7),'');
fprintf('   %-12s = %8.4f %s\n','C_Y',Y_eq_Ex4(8),'');
fprintf('   %-12s = %8.4f %s\n','C_l',Y_eq_Ex4(9),'');
fprintf('   %-12s = %8.4f %s\n','C_n',Y_eq_Ex4(10),'');
fprintf('   %-12s = %8.4f %s\n','rho',Y_eq_Ex4(11),'kg/m^3');
fprintf('   %-12s = %8.1f %s\n','qbar',Y_eq_Ex4(12),'N/m^2');
fprintf('\n');

%% Ex.5

X0_Ex5 = X_eq_Ex2;
U0_Ex5 = U_eq_Ex2;

dt_Ex5 = 0.010;
t_f_Ex5 = 30;
T_Ex5 = 0:dt_Ex5:t_f_Ex5;
[X_Ex5,Y_Ex5] = ode4xy(@simulate_flight,T_Ex5,X0_Ex5,U0_Ex5);
U_Ex5 = Y_Ex5(:,end-4:end);
Y_Ex5 = Y_Ex5(:,1:end-5);

plot_long
plot_latdir
plot_controls
plot_path
plot_outputs

%% Ex.6

x_eq_0_Ex6 = zeros(6,1);
x_eq_Ex6 = fsolve(@trim_function_Ex6,x_eq_0_Ex6,options,trim_par);
[~,X_eq_Ex6,U_eq_Ex6,Y_eq_Ex6] = trim_function_Ex6(x_eq_Ex6,trim_par);

nX_Ex6 = length(X_eq_Ex6);
nU_Ex6 = length(U_eq_Ex6);
nY_Ex6 = length(Y_eq_Ex6);

lin_output_Ex6 = struct('A',zeros(nX_Ex6,nX_Ex6),'B',zeros(nX_Ex6,nU_Ex6),'C',zeros(nY_Ex6,nX_Ex6),'D',zeros(nY_Ex6,nU_Ex6));

delta_val_Ex6 = 1e-5;

A_Ex6 = zeros(nX_Ex6,nX_Ex6);
C_Ex6 = zeros(nY_Ex6,nX_Ex6);
for j=1:nX_Ex6
    dX = zeros(nX_Ex6,1);
    dX(j) = delta_val_Ex6;
    [Xdot_plus, Y_plus]= dynamics_Ex6(0, X_eq_Ex6 + dX, U_eq_Ex6);
    [Xdot_minus, Y_minus]= dynamics_Ex6(0, X_eq_Ex6 - dX, U_eq_Ex6);
    A_Ex6(:, j) = (Xdot_plus - Xdot_minus)/(2*delta_val_Ex6);
    C_Ex6(:, j) = (Y_plus - Y_minus)/(2*delta_val_Ex6);
end

B_Ex6 = zeros(nX_Ex6,nU_Ex6);
D_Ex6 = zeros(nY_Ex6,nU_Ex6);
for j=1:nU_Ex6
    dU = zeros(nU_Ex6,1);
    dU(j) = delta_val_Ex6;
    [Xdot_plus, Y_plus]= dynamics_Ex6(0, X_eq_Ex6, U_eq_Ex6 + dU);
    [Xdot_minus, Y_minus]= dynamics_Ex6(0, X_eq_Ex6, U_eq_Ex6 - dU);
    B_Ex6(:, j) = (Xdot_plus - Xdot_minus)/(2*delta_val_Ex6);
    D_Ex6(:, j) = (Y_plus - Y_minus)/(2*delta_val_Ex6);
end

lin_output_Ex6.A = A_Ex6; lin_output_Ex6.B = B_Ex6; lin_output_Ex6.C = C_Ex6; lin_output_Ex6.D = D_Ex6;

sel_states_Ex6 = 1:12;
[eigvec_Ex6,eigval_Ex6] = eig(lin_output_Ex6.A(sel_states_Ex6,sel_states_Ex6));
eigval_Ex6 = diag(eigval_Ex6);
[eigval_Ex6,i_sort_Ex6] = sort(eigval_Ex6);
eigvec_Ex6 = eigvec_Ex6(:,i_sort_Ex6);

fprintf('\n----- EX.6 -----\n\n');
damp(eigval_Ex6)

[~,i_roll] = max(abs(eigvec_Ex3(9,:)));   
lam = eigval_Ex3(i_roll);
wn  = abs(lam); 
z   = -real(lam)/max(wn,eps); 
tau = -1/real(lam);
fprintf('EX.3 roll: λ=% .4f%+.4fi, ωn=%.3f, ζ=%.3f, τ=%.6fs\n', real(lam), imag(lam), wn, z, tau);

eig_roll = eigvec_Ex3(:,i_roll)/eigvec_Ex3(9,i_roll);
fprintf('EX.3 roll eigenvector (normalizado com p=1):\n');
fprintf('   beta   : % .6e\n', eig_roll(7));
fprintf('   phi    : % .6e\n', eig_roll(8));
fprintf('   p      : % .6e\n', eig_roll(9));
fprintf('   r      : % .6e\n', eig_roll(10));

[~,i_roll] = max(abs(eigvec_Ex6(9,:)));  
lam = eigval_Ex6(i_roll);
wn  = abs(lam); 
z   = -real(lam)/max(wn,eps); 
tau = -1/real(lam);
fprintf('EX.6 roll: λ=% .4f%+.4fi, ωn=%.3f, ζ=%.3f, τ=%.6fs\n', real(lam), imag(lam), wn, z, tau);

eig_roll6 = eigvec_Ex6(:,i_roll)/eigvec_Ex6(9,i_roll);
fprintf('EX.6 roll eigenvector (normalizado com p=1):\n');
fprintf('   beta   : % .6e\n', eig_roll6(7));
fprintf('   phi    : % .6e\n', eig_roll6(8));
fprintf('   p      : % .6e\n', eig_roll6(9));
fprintf('   r      : % .6e\n', eig_roll6(10));

%% Ex.7

nX_Ex7 = length(X_eq_Ex2);
nU_Ex7 = length(U_eq_Ex2);
nY_Ex7 = length(Y_eq_Ex2);

lin_output_Ex7 = struct('A',zeros(nX_Ex7,nX_Ex7),'B',zeros(nX_Ex7,nU_Ex7),'C',zeros(nY_Ex7,nX_Ex7),'D',zeros(nY_Ex7,nU_Ex7));

delta_val_Ex7 = 1e-5;

A_Ex7 = zeros(nX_Ex7,nX_Ex7);
C_Ex7 = zeros(nY_Ex7,nX_Ex7);
for j=1:nX_Ex7
    dX = zeros(nX_Ex7,1);
    dX(j) = delta_val_Ex7;
    [Xdot_plus, Y_plus]= dynamics(0, X_eq_Ex2 + dX, U_eq_Ex2);
    [Xdot_minus, Y_minus]= dynamics(0, X_eq_Ex2 - dX, U_eq_Ex2);
    A_Ex7(:, j) = (Xdot_plus - Xdot_minus)/(2*delta_val_Ex7);
    C_Ex7(:, j) = (Y_plus - Y_minus)/(2*delta_val_Ex7);
end

B_Ex7 = zeros(nX_Ex7,nU_Ex7);
D_Ex7 = zeros(nY_Ex7,nU_Ex7);
for j=1:nU_Ex7
    dU = zeros(nU_Ex7,1);
    dU(j) = delta_val_Ex7;
    [Xdot_plus, Y_plus]= dynamics(0, X_eq_Ex2, U_eq_Ex2 + dU);
    [Xdot_minus, Y_minus]= dynamics(0, X_eq_Ex2, U_eq_Ex2 - dU);
    B_Ex7(:, j) = (Xdot_plus - Xdot_minus)/(2*delta_val_Ex7);
    D_Ex7(:, j) = (Y_plus - Y_minus)/(2*delta_val_Ex7);
end

lin_output_Ex7.A = A_Ex7; lin_output_Ex7.B = B_Ex7; lin_output_Ex7.C = C_Ex7; lin_output_Ex7.D = D_Ex7;

sel_states_Ex7 = [7 8 9 10];
[eigvec_Ex7,eigval_Ex7] = eig(lin_output_Ex7.A(sel_states_Ex7, sel_states_Ex7));
eigval_Ex7 = diag(eigval_Ex7);
[eigval_Ex7,i_sort_Ex7] = sort(eigval_Ex7);
eigvec_Ex7 = eigvec_Ex7(:, i_sort_Ex7);

fprintf('\n----- EX.7 -----\n\n');
damp(eigval_Ex7)

Afull = lin_output_Ex7.A;

[VF,DF] = eig(Afull); lamF = diag(DF);
is_cpx = abs(imag(lamF)) > 1e-8;                     
scoreF = abs(VF(7,is_cpx)) + abs(VF(10,is_cpx));     
[~,kF] = max(scoreF);
lamF_dr = lamF(find(is_cpx,1,'first')-1 + kF);        
if imag(lamF_dr) < 0, lamF_dr = conj(lamF_dr); end
wnF = abs(lamF_dr); zF = -real(lamF_dr)/max(wnF,eps);

lamR_all = eig(lin_output_Ex7.A([7 8 9 10],[7 8 9 10]));
[~,kR] = max(abs(imag(lamR_all)));
lamR = lamR_all(kR); if imag(lamR)<0, lamR = conj(lamR); end
wnR = abs(lamR); zR = -real(lamR)/max(wnR,eps);
fprintf('FULL DR: ωn=%.4f ζ=%.4f | RED 4x4 DR: ωn=%.4f ζ=%.4f | Δωn=%.1f%% Δζ=%.1f%%\n', ...
    wnF, zF, wnR, zR, 100*(wnR-wnF)/wnF, 100*(zR-zF)/max(zF,eps));