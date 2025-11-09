function[f,X,U,Y] = trim_function(x,trim_par)

global aircraft

alpha_deg_eq  = x(1);
delta_t_eq    = x(2);
i_t_deg_eq    = x(3);

V_eq = trim_par.V;
h_eq = trim_par.h;

aircraft.i_t = i_t_deg_eq * (pi/180);

theta_deg_eq = alpha_deg_eq; 

X = [
     V_eq
     alpha_deg_eq
     0
     theta_deg_eq
     h_eq
     0
     0
     0
     0
     0
     0
     0
     ];

U = [
    delta_t_eq
    delta_t_eq
    0
    0
    0
    ];

[Xdot, Y] = dynamics(0, X, U);

f = [
    Xdot(1)
    Xdot(2)
    Xdot(3)
    ];
end