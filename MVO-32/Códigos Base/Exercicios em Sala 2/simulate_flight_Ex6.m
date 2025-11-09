function [Xdot,Y] = simulate_flight_Ex6(t,X,U,i_cond)

% Rudder doublet:
if t>=1 && t<3
    U(5) = U(5) - 5;
elseif t>=3 && t<5
    U(5) = U(5) + 5;
end

[Xdot,Y] = dynamics(t,X,U,i_cond);

Y = [
    Y
    U
    ];