clear
close all
clc

% MVO-32 ESTABILIDADE E CONTROLE DE AERONAVES
% Dupla: Thiago Akira Missato e Eduardo Pinto Ferreira.
% Código MAIN.m referente a Questão 03 da Lista 03.
% Elaborado por Eduardo Pinto Ferreira.


global g
global aircraft

g = 9.80665;


create_aircraft; 
load trim_output.mat; 

% linearização
i_cond = 1; 
nX = 12;
nU = 6;
nY = 12;
X_eq = trim_output(i_cond).X_eq;
U_eq = trim_output(i_cond).U_eq;
lin_output(i_cond) = struct('A',zeros(nX,nX),...
    'B',zeros(nX,nU),...
    'C',zeros(nY,nX),...
    'D',zeros(nY,nU));
delta_val = 1e-5; 

% jacobianas
A = zeros(nX,nX);
C = zeros(nY,nX);
for j=1:nX
    dX = zeros(nX,1);
    dX(j) = delta_val;
    [Xdot_plus,Y_plus] = dynamics(0,X_eq+dX,U_eq, i_cond);
    [Xdot_minus,Y_minus] = dynamics(0,X_eq-dX,U_eq, i_cond);
    A(:,j) = (Xdot_plus - Xdot_minus)/(2*dX(j));
    C(:,j) = (Y_plus - Y_minus)/(2*dX(j));
end
B = zeros(nX,nU);
D = zeros(nY,nU);
for j=1:nU
    dU = zeros(nU,1);
    dU(j) = delta_val;
    [Xdot_plus,Y_plus] = dynamics(0,X_eq,U_eq+dU, i_cond);
    [Xdot_minus,Y_minus] = dynamics(0,X_eq,U_eq-dU, i_cond);
    B(:,j) = (Xdot_plus - Xdot_minus)/(2*dU(j));
    D(:,j) = (Y_plus - Y_minus)/(2*dU(j));
end
lin_output(i_cond).A = A;
lin_output(i_cond).B = B;
lin_output(i_cond).C = C;
lin_output(i_cond).D = D;

save lin_output.mat lin_output

% Modos sem(x, y:
sel_states = [1 2 3 4 5 7 8 9 10 11]; % V, alpha, q, theta, h, beta, phi, p, r, psi
state_names = {'V'; 'alpha'; 'q'; 'theta'; 'h'; 'beta'; 'phi'; 'p'; 'r'; 'psi'};
A_reduced = A(sel_states, sel_states);
[eigvec, eigval_matrix] = eig(A_reduced);
eigval = diag(eigval_matrix);
fprintf('\n');
% Zeta e Tau
damp(eigval);

fprintf('\n Autovetores\n');
disp(state_names); 
fprintf('\n Autovalores\n');
disp(eigval'); 
fprintf('\nAutovetores\n');
disp(eigvec);
