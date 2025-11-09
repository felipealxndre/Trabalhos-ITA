function [f, X, U, Y] = trim_function_failure(x, trim_par)

alpha_deg = x(1);
delta1   = x(2);  
i_t_deg   = x(3);
delta_e_deg = 0;
delta_a_deg = x(5);
delta_r_deg = x(6);

theta_deg = alpha_deg;

X = [
    trim_par.V   % V [m/s]
    alpha_deg    % alpha [deg]
    0            % q [deg/s]
    theta_deg    % theta [deg]
    trim_par.h   % h [m]
    0            % x [m]
    0            % beta [deg]
    0            % phi [deg]
    0            % p [deg/s]
    0            % r [deg/s]
    0            % psi [deg]
    0            % y [m]
];

delta2 = 0.15;

U = [
    delta1
    delta2
    i_t_deg
    delta_e_deg
    delta_a_deg
    delta_r_deg
];

[Xdot, Y] = dynamics(0, X, U);

V_eq = trim_par.V;
gamma_deg_eq = trim_par.gamma_deg;

hdot_eq = V_eq*sind(gamma_deg_eq);
xdot_eq = V_eq*cosd(gamma_deg_eq);

f = [
    Xdot(1:3)
    Xdot(4) - trim_par.thetadot_deg_s
    Xdot(5) - hdot_eq
    Xdot(6) - xdot_eq
    ];

end
