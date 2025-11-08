function aircraft = create_aircraft()
% This script creates the structure variable 'aircraft'

% DO NOT change the fields' names!

%--------------------------------------------------------------------------
% Geometric properties:
aircraft.S = 81.25;              % Reference wing area [m^2] #ok!
aircraft.c = 2.53;               % Wing mean aerodynamic chord [m] #ok!
aircraft.b = 35.99;              % Wing span [m] #ok!

%--------------------------------------------------------------------------
% Inertia properties:
aircraft.m = 41347.03;           % Mass [kg] #ok!
aircraft.Ixx = 474148.04;        % I_{xx} [kg.m^2] #ok!
aircraft.Iyy = 2527894.43;       % I_{yy} [kg.m^2] #ok!
aircraft.Izz = 2948507.61;       % I_{zz} [kg.m^2] #ok!
aircraft.Ixz = 19981.74;         % I_{xz} [kg.m^2] #ok!

%--------------------------------------------------------------------------
% Propulsive model properties:
aircraft.iota_1_deg = 2.00;      % Left engine incidence [deg] #ok!
aircraft.tau_1_deg = 1.50;       % Left engine toe-in angle [deg] #ok!
aircraft.x_1 = 13.90;            % Left engine thrust action point x_b coordinate [m] #ok!
aircraft.y_1 = -5.00;            % Left engine thrust action point y_b coordinate [m] #ok!
aircraft.z_1 = 1.00;             % Left engine thrust action point z_b coordinate [m] #ok!

aircraft.iota_2_deg = 2.00;      % Right engine incidence [deg] #ok!
aircraft.tau_2_deg = -1.50;       % Right engine toe-in angle [deg] #ok!
aircraft.x_2 = 13.90;            % Right engine thrust action point x_b coordinate [m] #ok!
aircraft.y_2 = 5.00;             % Right engine thrust action point y_b coordinate [m] #ok!
aircraft.z_2 = 1.00;             % Right engine thrust action point z_b coordinate [m] #ok!

aircraft.Tmax = 116445.95;       % Maximum sea-level thrust for a single engine [N] #ok!
aircraft.n_rho = 0.8;            % Air density exponent in the thrust model [-] #ok!

%--------------------------------------------------------------------------
% Aerodynamic data bank:

% Important notes:
% -> Aerodynamic moments calculated about the nominal CG, in the
% directions of the body axes
% -> Angular velocity components refer to rotational motion about the
% nominal CG
% -> All angles used in the aerodynamic databank will need to be converted 
% to radians
% -> All angular velocities used in the aerodynamic databank will need to 
% be converted to radians per second and then made dimensionless by the
% appropriate multiplication by c/(2*V) or b/(2*V)

% Lift coefficient:
aircraft.CL_0 = -0.0407;                             % [-] #ok!
aircraft.CL_alpha = 8.44;                            % [1/rad] #ok!
aircraft.CL_q = 19.4;                                % [1/rad] #ok!
aircraft.CL_i_t = 1.33;                              % [1/rad] #ok!
aircraft.CL_delta_e = 0.911;                         % [1/rad] #ok!

% Pitching moment coefficient:
aircraft.Cm_0 = -0.0596;                             % [-] #ok!
aircraft.Cm_alpha = -3.06;                           % [1/rad] #ok!
aircraft.Cm_q = -110.0;                              % [1/rad] #ok!
aircraft.Cm_i_t = -7.82;                             % [1/rad] #ok!
aircraft.Cm_delta_e = -5.60;                         % [1/rad] #ok!

% Drag coefficient:
aircraft.CD_0 = 0.0190;                              % [-] #ok!
aircraft.CD_alpha = 0.0;                             % [1/rad] #ok!
aircraft.CD_alpha2 = 0.0005;                         % [1/rad^2] #ok!
aircraft.CD_q = 0.516;                               % [1/rad] #ok!
aircraft.CD_i_t_0 = -0.0256;                         % [1/rad] #ok!
aircraft.CD_i_t_alpha = 0.3914;                      % [1/rad^2] #ok!
aircraft.CD_delta_e_0 = -0.0175;                     % [1/rad] #ok!
aircraft.CD_delta_e_alpha = 0.2677;                  % [1/rad^2] #ok!
aircraft.CD_beta2 = -0.55282;                         % [1/rad^2] #ok!
aircraft.CD_p2 = -1.244;                             % [1/rad^2] #ok!
aircraft.CD_r2 = -0.432;                             % [1/rad^2] #ok!
aircraft.CD_delta_a2 = 0.2284;                       % [1/rad^2] #ok!
aircraft.CD_delta_r2 = 0.1234;                       % [1/rad^2] #ok!

% Side force coefficient:
aircraft.CY_beta = 0.857;                            % [1/rad] #ok!
aircraft.CY_p_0 = 0.226;                             % [1/rad] #ok!
aircraft.CY_p_alpha = -4.9023;                       % [1/rad^2] #ok!
aircraft.CY_r = -0.0799;                             % [1/rad] #ok!
aircraft.CY_delta_a = 0.1197;                           % [1/rad] #ok! (AVL)
aircraft.CY_delta_r = -0.192;                        % [1/rad] #ok!

% Rolling moment coefficient:
aircraft.Cl_beta_0 = -2.00e-6;                       % [1/rad] #ok!
aircraft.Cl_beta_alpha = 3.0e-5;                     % [1/rad^2] #ok!
aircraft.Cl_p = -0.716;                              % [1/rad] #ok!
aircraft.Cl_r_0 = 0.0618;                            % [1/rad] #ok!
aircraft.Cl_r_alpha = 0.9504;                        % [1/rad^2] #ok!
aircraft.Cl_delta_a = -0.550;                        % [1/rad] #ok!
aircraft.Cl_delta_r = 0.0215;                        % [1/rad] #ok!

% Yawing moment coefficient:
aircraft.Cn_beta = 0.00445;                         % [1/rad] #ok! Sinal contrario do avl
aircraft.Cn_p_0 = 0.0289;                            % [1/rad] #ok!
aircraft.Cn_p_alpha = -1.5474;                       % [1/rad^2] #ok!
aircraft.Cn_r = -0.223;                              % [1/rad] #ok!
aircraft.Cn_delta_a_0 = -0.012;                      % [1/rad] #ok!
aircraft.Cn_delta_a_alpha = 0.1881;                  % [1/rad^2] #ok!
aircraft.Cn_delta_r = -0.0939;                       % [1/rad] #ok!

end