function [f,X,U,Y] = trim_function(x, trim_par, i_cond)
% Modificado para Lista 3 / Ex 2:
% 1. Aceita i_cond para compatibilidade
% 2. Implementa hipótese de empuxo simétrico: delta_1 = delta_2 = x(9)
% 3. Implementa restrição de profundor nulo: delta_e = 0

% --- Restrições para voo simétrico e do Ex 2 ---
beta_deg_eq = 0; % Voo simétrico
delta_e_deg_eq = 0; % Restrição do Ex 2 (equilíbrio com i_t)

% --- Vetor de Estado X (12x1) ---
% Incógnitas são x(1) a x(8)
X = [
     x(1) % V
     x(2) % alpha
     x(3) % q
     x(4) % theta
     trim_par.h % h (fixo)
     0 % x (fixo)
     beta_deg_eq % beta (fixo)
     x(5) % phi
     x(6) % p
     x(7) % r
     x(8) % psi
     0 % y (fixo)
    ];

% --- Vetor de Controle U (6x1) ---
% Incógnitas são x(9) a x(12)
U = [
    x(9) % U(1) = delta_1 (incógnita "throttle simétrico")
    x(9) % U(2) = delta_2 (Hipótese: delta_2 = delta_1)
    x(10) % U(3) = i_t (incógnita)
    delta_e_deg_eq % U(4) = delta_e (Restrição: delta_e = 0)
    x(11) % U(5) = delta_a (incógnita)
    x(12) % U(6) = delta_r (incógnita)
    ];

% --- Chama a dinâmica ---
% CORREÇÃO: Passando i_cond para dynamics
[Xdot,Y] = dynamics(0,X,U, i_cond);

% --- Alvos da condição de equilíbrio (para voo reto nivelado) ---
V_eq = trim_par.V;
gamma_deg_eq = trim_par.gamma_deg;
ydot_eq = 0;
hdot_eq = V_eq*sind(gamma_deg_eq); % hdot = 0 para gamma=0
xdot_eq = V_eq*cosd(gamma_deg_eq); % xdot = V para gamma=0

% --- Vetor de Resíduos f (12x1) ---
% O fsolve tentará zerar este vetor
f = [
    Xdot(1) % Vdot = 0
    Xdot(2) % alphadot = 0
    Xdot(3) % qdot = 0
    Xdot(4) - trim_par.thetadot_deg_s % thetadot = 0
    Xdot(5) - hdot_eq % hdot = 0
    Xdot(6) - xdot_eq % xdot = V
    Xdot(7) % betadot = 0
    Xdot(8) % phidot = 0
    Xdot(9) % pdot = 0
    Xdot(10) % rdot = 0
    Xdot(11) - trim_par.psidot_deg_s % psidot = 0
    Xdot(12) - ydot_eq % ydot = 0
    ];

end