function [f, X, U, Y] = trim_function_Ex6(x, trim_par)

alpha_deg   = x(1);
delta_c     = x(2);      % δ1 = δ2 = δc (simetria em cruzeiro)
it_deg      = x(3);

delta_e_deg = 0;         
delta_a_deg = 0;        
delta_r_deg = 0;      

V         = trim_par.V;
h         = trim_par.h;
theta_deg = alpha_deg;
q_deg_s   = 0;
beta_deg  = 0;
phi_deg   = 0;
p_deg_s   = 0;
r_deg_s   = 0;
psi_deg   = 0;
x_pos     = 0;
y_pos     = 0;

X = [ V; alpha_deg; q_deg_s; theta_deg; h; x_pos; beta_deg; phi_deg; ...
      p_deg_s; r_deg_s; psi_deg; y_pos];

delta1 = delta_c;
delta2 = delta_c;
U = [ delta1; delta2; it_deg; delta_e_deg; delta_a_deg; delta_r_deg ];

[Xdot, Y] = dynamics_Ex6(0, X, U); 

V_eq         = trim_par.V;
gamma_deg_eq = 0;
thetadot_sp  = 0;
hdot_eq      = V_eq * sind(gamma_deg_eq);  
xdot_eq      = V_eq * cosd(gamma_deg_eq);  

f = [ Xdot(1:3);
      Xdot(4)  - thetadot_sp;
      Xdot(5)  - hdot_eq;
      Xdot(6)  - xdot_eq ];
end
