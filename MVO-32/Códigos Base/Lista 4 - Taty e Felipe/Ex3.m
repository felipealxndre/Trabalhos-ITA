clear
clc

[A,B,C,D,X_eq,U_eq,Y_eq] = gera_modelo_L3();

tau_throttle = 0.1;

A_base = A;
B_base = B;
C_base = C;
D_base = D;

[nx,nu] = size(B_base);
ny = size(C_base,1);

B_delta = B_base(:,1:2);
B_rest  = B_base(:,3:end);

A_throttle = zeros(nx+2,nx+2);
A_throttle(1:nx,1:nx) = A_base;
A_throttle(1:nx,nx+1:nx+2) = B_delta;
A_throttle(nx+1,nx+1) = -1/tau_throttle;
A_throttle(nx+2,nx+2) = -1/tau_throttle;

B_throttle = zeros(nx+2,nu);
B_throttle(1:nx,3:nu) = B_rest;
B_throttle(nx+1,1) = 1/tau_throttle;
B_throttle(nx+2,2) = 1/tau_throttle;

C_throttle = [C_base zeros(ny,2)];
D_throttle = D_base;

fprintf('\nAutovalores SEM dinâmicas dos motores:\n\n');
damp(A_base)

fprintf('\nAutovalores COM dinâmicas dos motores:\n\n');
damp(A_throttle)

lambda_com = eig(A_throttle);
polos_reais = lambda_com(abs(imag(lambda_com)) < 1e-6);
polos_reais = sort(real(polos_reais));

fprintf('\nPolos reais do sistema com dinâmica de throttle:\n');
for k = 1:length(polos_reais)
    fprintf('   % .4f\n', polos_reais(k));
end
