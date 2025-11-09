function [Xdot,Y] = simulate_flight_Ex5(t,X,Ucom,i_cond)

tau_e = 0.050;

delta_e_deg = X(7);

% Elevator doublet:
if t>=1 && t<3
    Ucom(3) = Ucom(3) + 3;
elseif t>=3 && t<5
    Ucom(3) = Ucom(3) - 3;
end

u_e_deg = Ucom(3);

delta_e_dot_deg_s = -1/tau_e*delta_e_deg + 1/tau_e*u_e_deg;

U = [
    Ucom(1:2)
    delta_e_deg
    ];


[Xdot,Y] = long_dynamics(t,X(1:6),U,i_cond);

Xdot = [
    Xdot
    delta_e_dot_deg_s
    ];

Y = [
    Y
    U
    u_e_deg
    ];