function[f,X,U,Y] = trim_function(x,trim_par)

%x = [V; alpha_deg; q_deg_s; theta_deg; throttle; i_t_deg];

delta_e_deg_eq = 0;

X = [
     x(1:4)  %V; alpha_deg; q_deg_s; theta_deg;
     trim_par.h
     0
     ];

U = [
    x(5:6) %throttle; i_t_deg
    delta_e_deg_eq
    ];

[Xdot, Y] = dynamics(0, X, U, i_cond);


V_eq = trim_par.V;
gamma_deg_eq = trim_par.gamma_deg;

hdot_eq = V_eq*sind(gamma_deg_eq);
xdot_eq = V_eq*cosd(gamma_deg_eq);

f = [
    Xdot(1:3) %Vdot; alphadot_deg_s; qdot_deg_s_2
    Xdot(4) - trim_par.thetadot_deg_s
    Xdot(5) - hdot_eq
    Xdot(6) - xdot_eq
    ];
