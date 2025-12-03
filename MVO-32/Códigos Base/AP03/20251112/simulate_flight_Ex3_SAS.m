function [Xdot,Y] = simulate_flight_Ex3_SAS(t,X,U_pilot,i_cond)

global aircraft

% Elevator doublet:
if t>=1 && t<3
    U_pilot(3) = U_pilot(3) + 3;
elseif t>=3 && t<5
    U_pilot(3) = U_pilot(3) - 3;
end


k_alpha = aircraft.SAS.k_alpha;
k_q = aircraft.SAS.k_q;
alpha_deg_eq = aircraft.SAS.alpha_deg_eq;
q_deg_s_eq = aircraft.SAS.q_deg_s_eq;

alpha_deg = X(2);
q_deg_s = X(3);

Delta_alpha_deg = alpha_deg - alpha_deg_eq;
Delta_q_deg_s = q_deg_s - q_deg_s_eq;

Delta_delta_e_deg_SAS = -k_alpha*Delta_alpha_deg - k_q*Delta_q_deg_s;

U = U_pilot;
U(3) = U(3) + Delta_delta_e_deg_SAS;


[Xdot,Y] = dynamics(t,X,U,i_cond);

Y = [
    Y
    U
    Delta_delta_e_deg_SAS
    ];