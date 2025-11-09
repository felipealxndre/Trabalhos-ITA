function [Xdot,Y] = simulate_flight_Ex5(t,X,U,i_cond)

% Elevator doublet:
if t>=1 && t<3
    U(3) = U(3) + 3;
elseif t>=3 && t<5
    U(3) = U(3) - 3;
end

[Xdot,Y] = long_dynamics(t,X,U,i_cond);

Y = [
    Y
    U
    ];