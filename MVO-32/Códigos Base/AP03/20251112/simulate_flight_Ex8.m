function [Xdot,Y] = simulate_flight_Ex8(t,X,U,i_cond)

[Xdot,Y] = dynamics_Ex8(t,X,U,i_cond);

Y = [
    Y
    U
    ];