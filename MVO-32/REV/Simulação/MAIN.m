clear
close all
clc

global g
global aircraft

g = 9.80665;

aircraft = struct('S', 40, 'c', 2.04, 'b', 20.8, ...
    'm', 10500, ...
    'Ixx', 1.03e5, 'Iyy', 1.58e5, 'Izz', 2.40e5, 'Ixz', 1.25e4, ...
    'i_p_deg', -1.5, ...
    'x_p', 1.50, 'z_p', -0.55, ...
    'Tmax', 27350, ...
    'n_rho', 0.8, ...
    'n_V', -1, ...
    'V_ref', 120.49);

psidot_deg_s_eq = 0;

trim_par(1) = struct('V', 120.49, 'h', 5486.4, 'gamma_deg', 0, ...
    'thetadot_deg_s', 0, 'psidot_deg_s', psidot_deg_s_eq);

options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);

trim_output(1) = struct('X_eq',zeros(12,1),'U_eq',zeros(4,1),'Y_eq',zeros(11,1));

for i_cond = 1:1
    x_eq_0 = zeros(12,1);
    x_eq_0(1) = trim_par(i_cond).V;
    x_eq_0(2) = 2;
    x_eq_0(4) = 2;
    x_eq_0(9) = 0.4;
    
    x_eq = fsolve(@trim_function, x_eq_0, options, trim_par, i_cond);
    
    [~,X_eq,U_eq,Y_eq] = trim_function(x_eq, trim_par, i_cond);
    
    trim_output(i_cond).X_eq = X_eq;
    trim_output(i_cond).U_eq = U_eq;
    trim_output(i_cond).Y_eq = Y_eq;

    fprintf('----- CONDIÇÃO DE VOO %d -----\n\n',i_cond);
    fprintf('   %-10s = %10.2f %-4s\n','V',X_eq(1),'m/s');
    fprintf('   %-10s = %10.4f %-4s\n','alpha',X_eq(2),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','theta',X_eq(4),'deg');
    fprintf('   %-10s = %10.1f %-4s\n','h',X_eq(5),'m');
    fprintf('\n');
    fprintf('   %-10s = %10.4f %-4s\n','beta',X_eq(7),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','phi',X_eq(8),'deg');
    fprintf('\n');
    fprintf('   %-10s = %10.2f %-4s\n','throttle',U_eq(1)*100,'%');
    fprintf('   %-10s = %10.4f %-4s\n','delta_e',U_eq(2),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','delta_a',U_eq(3),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','delta_r',U_eq(4),'deg');
    fprintf('\n');
end

save trim_output.mat trim_output

i_cond = 1;
X_trim = trim_output(i_cond).X_eq;
U_trim = trim_output(i_cond).U_eq;

t_final = 10;
dt = 0.01;
t_span = 0:dt:t_final;

[T_sim, X_sim] = ode45(@(t,x) dynamics_with_input(t, x, U_trim, i_cond), t_span, X_trim);

figure('Name', 'Resposta do Período Curto');
subplot(3,1,1);
plot(T_sim, X_sim(:,2));
grid on; ylabel('\alpha [deg]'); title('Ângulo de Ataque');
subplot(3,1,2);
plot(T_sim, X_sim(:,3));
grid on; ylabel('q [deg/s]'); title('Taxa de Arfagem');
subplot(3,1,3);
plot(T_sim, X_sim(:,1));
grid on; ylabel('V [m/s]'); xlabel('Tempo [s]'); title('Velocidade');

function [dx] = dynamics_with_input(t, x, U_trim, i_cond)
    U_current = U_trim;
    amplitude = 2;
    t_start = 1.0;
    duration = 1.0;
    if t >= t_start && t < (t_start + duration/2)
        U_current(2) = U_trim(2) - amplitude;
    elseif t >= (t_start + duration/2) && t < (t_start + duration)
        U_current(2) = U_trim(2) + amplitude;
    end
    [dx, ~] = dynamics(t, x, U_current, i_cond);
end

i_cond = 1;
X_eq = trim_output(i_cond).X_eq;
U_eq = trim_output(i_cond).U_eq;

nX = length(X_eq);
perturbacao = 1e-5;
A = zeros(nX, nX);

for k = 1:nX
    X_plus = X_eq; X_plus(k) = X_plus(k) + perturbacao;
    X_minus = X_eq; X_minus(k) = X_minus(k) - perturbacao;
    [dX_plus, ~] = dynamics(0, X_plus, U_eq, i_cond);
    [dX_minus, ~] = dynamics(0, X_minus, U_eq, i_cond);
    A(:, k) = (dX_plus - dX_minus) / (2 * perturbacao);
end

idx_long = [1, 2, 3, 4];
A_long = A(idx_long, idx_long);

autovalores = eig(A_long);

fprintf('\n--- AUTOVALORES LONGITUDINAIS ---\n');
disp(autovalores);

fprintf('\n--- CARACTERÍSTICAS DOS MODOS ---\n');
damp(autovalores);
