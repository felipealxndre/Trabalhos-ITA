function [Xdot,Y] = simulate_flight_Ex5(t,X,U,i_cond)

% Doublet de empuxo diferencial nos motores
if t>=1 && t<3.5
    U(1) = U(1) - 0.1;
    U(2) = U(2) + 0.1;
elseif t>=3.5 && t<6
    U(1) = U(1) + 0.1;
    U(2) = U(2) - 0.1;
end

[Xdot,Y] = dynamics(t,X,U,i_cond);

Y = [
    Y
    U
    ];