function [Xdot,Y] = simulate_flight(t,X,U)

delta1 = U(1);
delta2 = U(2);

% Rudder doublet:
if t>=1 && t<3.5
    U(1) = delta1 - 0.1;
    U(2) = delta2 + 0.1;
elseif t>=3.5 && t<6
    U(1) = delta1 + 0.1;
    U(2) = delta2 - 0.1;
end

[Xdot,Y] = dynamics(t,X,U);

Y = [
    Y
    U
    ];