function [Xdot,Y] = simulate_flight(t,X,U,i_cond)

[Xdot,Y] = long_dynamics(t,X,U,i_cond);

Y = [
    Y
    U
    ];