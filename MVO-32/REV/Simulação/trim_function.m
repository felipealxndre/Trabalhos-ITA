function [f, X, U, Y] = trim_function(x, trim_par, i_cond)

% Definição do vetor de incógnitas x (12 variáveis):
% x(1:8): [V; alpha; q; theta; phi; p; r; psi]
% x(9): Throttle
% x(10): Delta_e (Profundor) - Substitui o antigo i_t
% x(11): Delta_a
% x(12): Delta_r

% Dornier 328 não tem estabilizador móvel (trimável)
beta_deg_eq = 0;

X = [
    x(1:4)
    trim_par(i_cond).h
    0
    beta_deg_eq
    x(5:8)
    0
    ];

% Montagem do Vetor de Controle U
U = [
    x(9)      % Throttle
    x(10)     % delta_e (Profundor)
    x(11)     % delta_a
    x(12)     % delta_r
    ];

[Xdot, Y] = dynamics(0, X, U, i_cond);

V_eq = trim_par(i_cond).V;
gamma_deg_eq = trim_par(i_cond).gamma_deg;
ydot_deg_eq = 0;
hdot_eq = V_eq*sind(gamma_deg_eq);
xdot_eq = V_eq*cosd(gamma_deg_eq);


f = [
    Xdot(1:3)
    Xdot(4) - trim_par(i_cond).thetadot_deg_s
    Xdot(5) - hdot_eq
    Xdot(6) - xdot_eq
    Xdot(7:10)
    Xdot(11) - trim_par(i_cond).psidot_deg_s
    Xdot(12) - ydot_deg_eq
    ];
end