function [f, X, U, Y] = trim_function_Ex4(x, trim_par)
%
% Esta função calcula o equilíbrio 6-DOF para a condição de falha do Ex. 4.
% Ela é modelada na estratégia do código de referência, resolvendo
% um sistema 12x12.
%
% Assume-se equilíbrio com beta=0 (sem derrapagem), mas phi (inclinação)
% é livre para compensar as forças laterais.

<<<<<<< HEAD:MVO-32/Códigos Base/trim_function_Ex4.m
% 1. Definir valores fixos para esta trimagem
beta_deg_eq    = 0;       % Assunção de voo com zero derrapagem
delta_e_deg_eq = 0;       % Assunção de trimagem com i_t (do Ex. 2)
delta_2_eq     = 0.15;    % Condição de falha do Ex. 4
=======
alpha_deg    = x(1);
delta1 = min(x(2), 1);
it_deg       = x(3);  
phi_deg      = x(4);
delta_a_deg  = x(5);
delta_r_deg  = x(6);
>>>>>>> e73cdf8fed4e7dfb7bc6ddcb53851918cda23541:MVO-32/Códigos Base/Lista 3 - Taty e Felipe/trim_function_Ex4.m

% 2. Desempacotar as 12 incógnitas do vetor 'x'
%    x = [V, alpha, q, theta, phi, p, r, psi, delta_1, i_t, delta_a, delta_r]
V_in         = x(1);
alpha_deg_in = x(2);
q_deg_s_in   = x(3);
theta_deg_in = x(4);
phi_deg_in   = x(5);  % Incógnita (chave para o equilíbrio)
p_deg_s_in   = x(6);
r_deg_s_in   = x(7);
psi_deg_in   = x(8);
delta_1_in   = min(x(9), 1); % Saturação em 100%
i_t_deg_in   = x(10);
delta_a_in   = x(11);
delta_r_in   = x(12);

<<<<<<< HEAD:MVO-32/Códigos Base/trim_function_Ex4.m
% 3. Montar o vetor de estado X (12 estados)
X = [
    V_in;
    alpha_deg_in;
    q_deg_s_in;
    theta_deg_in;
    trim_par.h;     % h (fixo do trim_par)
    0;              % x_pos (irrelevante para o equilíbrio)
    beta_deg_eq;    % beta (fixo em 0)
    phi_deg_in;     % phi (incógnita)
    p_deg_s_in;
    r_deg_s_in;
    psi_deg_in;
    0               % y_pos (irrelevante para o equilíbrio)
];
=======
V            = trim_par.V;
theta_deg    = alpha_deg;
q_deg_s      = 0;
h            = trim_par.h;
beta_deg     = 0; 
p_deg_s      = 0;
r_deg_s      = 0;
psi_deg      = 0;
x_pos        = 0;
y_pos        = 0;
>>>>>>> e73cdf8fed4e7dfb7bc6ddcb53851918cda23541:MVO-32/Códigos Base/Lista 3 - Taty e Felipe/trim_function_Ex4.m

% 4. Montar o vetor de controle U (6 controles)
U = [
    delta_1_in;
    delta_2_eq;       % Fixo (falha)
    i_t_deg_in;
    delta_e_deg_eq;   % Fixo (suposição)
    delta_a_in;
    delta_r_in
];

% 5. Rodar a dinâmica
[Xdot, Y] = dynamics(0, X, U);

% 6. Definir o "Alvo" (condições de equilíbrio)
%    O alvo é voo em equilíbrio (sem acelerações), nivelado e reto.
V_eq         = V_in; % O 'V_eq' do trim_par era só um chute
gamma_deg_eq = 0;

% O vetor de alvo (o que Xdot DEVERIA ser)
Xdot_eq = zeros(12,1); % Quase tudo é zero

% Condições cinemáticas de voo reto nivelado
% (O código de referência faz isso de forma mais simples,
% subtraindo os alvos no vetor f)

% 7. Montar o vetor de erro 'f' (12 equações)
%    f = (O que está acontecendo) - (O que queremos)
%    Esta é a forma como o código de referência faz,
%    garantindo que todas as 12 derivadas sejam zeradas.

hdot_eq = V_eq * sind(gamma_deg_eq);
xdot_eq = V_eq * cosd(gamma_deg_eq);
ydot_eq = 0;
thetadot_eq = 0; % Para voo reto
psidot_eq = 0; % Para voo reto

f = [
    Xdot(1:3); % Vdot, alphadot, qdot (queremos 0)
    Xdot(4) - thetadot_eq; % thetadot
    Xdot(5) - hdot_eq;     % hdot
    Xdot(6) - xdot_eq;     % xdot
    Xdot(7:10);% betadot, phidot, pdot, rdot (queremos 0)
    Xdot(11) - psidot_eq;  % psidot
    Xdot(12) - ydot_eq     % ydot
];

end