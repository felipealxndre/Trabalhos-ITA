function [f, X, U, Y] = trim_function(x, trim_par)
% x = [ V; alpha_deg; q_deg_s; theta_deg; delta_c; i_t_deg ];
% Hipótese no equilíbrio (Ex. 2 L3):
%   - Voar reto nivelado: γ = 0, ψdot = 0
%   - Simetria propulsiva: δ1 = δ2 = δc
%   - δe = 0, δa = 0, δr = 0
%   - Sem acoplamentos laterais no trim: β = 0, ϕ = 0, p = 0, r = 0

% Desempacota incógnitas
V              = x(1);
alpha_deg      = x(2);
q_deg_s        = x(3);
theta_deg      = x(4);
delta_c        = x(5);   % throttle comum (δ1 = δ2 = δc)
i_t_deg        = x(6);


% Monta estados completos (ver Eq. (20) da Lista 3)
% X = [ V α q θ h x β ϕ p r ψ y ]^T  (ângulos em deg; vel ang em deg/s)
X = [
    V;
    alpha_deg;
    q_deg_s;
    theta_deg;
    trim_par.h;
    0;           % x (posição ao longo)
    0;           % beta
    0;           % phi
    0;           % p
    0;           % r
    0;           % psi
    0            % y (lateral)
];

% Monta controles (ver Eq. (21) da Lista 3)
% U = [ δ1 δ2 it δe δa δr ]^T
delta1 = delta_c;
delta2 = delta_c;
delta_e_deg = 0;
delta_a_deg = 0;
delta_r_deg = 0;

U = [
    delta1;
    delta2;
    i_t_deg;
    delta_e_deg;
    delta_a_deg;
    delta_r_deg
];

% Dinâmica completa do avião (deve refletir X, U e Y conforme a Lista 3)
% Y = [ γ T1 T2 Mach CD CL Cm CY Cl Cn ρ qbar ]^T   (Eq. (22))
[Xdot, Y] = dynamics(0, X, U);

% Condição-alvo: voo reto nivelado (γ = 0) e ψdot = 0
V_eq = trim_par.V;
hdot_eq = 0;           % γ = 0
xdot_eq = V_eq;        % ẋ = V (reta; ψdot = 0)
thetadot_eq = trim_par.thetadot_deg_s; % 0 deg/s

f = [
    Xdot(1);                     % Vdot = 0
    Xdot(2);                     % alphadot = 0
    Xdot(3);                     % qdot = 0
    Xdot(4) - thetadot_eq;       % thetadot = 0
    Xdot(5) - hdot_eq;           % hdot = 0
    Xdot(6) - xdot_eq            % xdot = V_eq
];
end
