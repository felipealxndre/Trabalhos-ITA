function [f, X,U,Y] = trim_function(x,trim_par,i_cond)

% x = [V; alpha_deg; q_deg_s; theta_deg; 
%       phi_deg; p_deg_s; r_deg_s; psi_deg;
%       throttle; i_t_deg; delta_a_deg; delta_r_deg] 

delta_e_deg_eq = 0;
beta_deg_eq = 0;

X = [
    x(1:4)
    trim_par(i_cond).h
    0
    beta_deg_eq
    x(5:8)
    0
    ];
U = [
    x(9:10)
    delta_e_deg_eq
    x(11:12)
    ];

[Xdot, Y] = dynamics(0,X,U,i_cond);

% Chapter 9, slide 18
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
    Xdot(11) - trim_par(i_cond ).psidot_deg_s
    Xdot(12) - ydot_deg_eq
    ];