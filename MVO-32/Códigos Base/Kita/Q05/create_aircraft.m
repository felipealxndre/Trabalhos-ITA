%   Posicao aft do CG 15.37
x_CG=-15.37;
z_CG=1.75;

% This script creates the structure variable 'aircraft'

% DO NOT change the fields' names!

%--------------------------------------------------------------------------
% Geometric properties:
aircraft.S = 99.02;              % Reference wing area [m^2]
aircraft.c = 3.574;              % Wing mean aerodynamic chord [m]
aircraft.b = 31.80;              % Wing span [m]

%--------------------------------------------------------------------------
% Inertia properties:
aircraft.m = 46192.12;              % Mass [kg]
aircraft.Ixx = 379542.17;            % I_{xx} [kg.m^2]
aircraft.Iyy = 2788136.18;            % I_{yy} [kg.m^2]
aircraft.Izz = 3073815.45;            % I_{zz} [kg.m^2]
aircraft.Ixz = 0;            % I_{xz} [kg.m^2]

%--------------------------------------------------------------------------
% Propulsive model properties:
aircraft.iota_1_deg = 2.0;     % Left engine incidence [deg]
aircraft.tau_1_deg = 1.5;      % Left engine toe-in angle [deg]
aircraft.x_1 = -11.164-x_CG;            % Left engine thrust action point x_b coordinate [m]
aircraft.y_1 = -4.2;            % Left engine thrust action point y_b coordinate [m]
aircraft.z_1 = 1.8715-z_CG;            % Left engine thrust action point z_b coordinate [m]

aircraft.iota_2_deg = 2.0;     % Right engine incidence [deg]
aircraft.tau_2_deg = -1.5;      % Right engine toe-in angle [deg]
aircraft.x_2 = -11.164-x_CG;            % Right engine thrust action point x_b coordinate [m]
aircraft.y_2 = 4.2;            % Right engine thrust action point y_b coordinate [m]
aircraft.z_2 = 1.8715-z_CG;            % Right engine thrust action point z_b coordinate [m]

% General Eletric CF34-8E Turbofans
aircraft.Tmax = 62300;           % Maximum sea-level thrust for a single engine [N]
aircraft.n_rho = 0.8;          % Air density exponent in the thrust model [-]

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
aircraft.CL_0 = 0.26222;                               % [-]
aircraft.CL_alpha = 7.140901;                           % [1/rad]
aircraft.CL_q = 16.757652;                               % [1/rad]
aircraft.CL_i_t = 0.93021;                             % [1/rad]
aircraft.CL_delta_e = 0.65031;                         % [1/rad]

% Pitching moment coefficient:
aircraft.Cm_0 = 0.07834;                               % [-]
aircraft.Cm_alpha = -3.247969;                           % [1/rad]
aircraft.Cm_q = -39.103592;                               % [1/rad]
aircraft.Cm_i_t = -3.9064;                             % [1/rad]
aircraft.Cm_delta_e = -2.8264;                         % [1/rad]

% Drag coefficient:
aircraft.CD_0 = 0.02172;                               % [-]
aircraft.CD_alpha = 0.17854;                           % [1/rad]
aircraft.CD_alpha2 = 1.61303;                          % [1/rad^2]
aircraft.CD_q = -0.186768;                               % [1/rad]
aircraft.CD_i_t_0 = 0.00905;                           % [1/rad]
aircraft.CD_i_t_alpha = 1.1335;                       % [1/rad^2]
aircraft.CD_delta_e_0 = 0.00613;                       % [1/rad]
aircraft.CD_delta_e_alpha = 0.8059;                   % [1/rad^2]
aircraft.CD_beta2 = -0.328;                           % [1/rad^2]
aircraft.CD_p2 = -0.828;                              % [1/rad^2]
aircraft.CD_r2 = -0.248;                              % [1/rad^2]
aircraft.CD_delta_a2 = 0.121;                        % [1/rad^2]
aircraft.CD_delta_r2 = 0.070;                        % [1/rad^2]

% Side force coefficient:
aircraft.CY_beta = 0.610322;                            % [1/rad]
aircraft.CY_p_0 = 0.122959;                             % [1/rad]
aircraft.CY_p_alpha = 13.141;                         % [1/rad^2]
aircraft.CY_r = -0.483905;                               % [1/rad]
aircraft.CY_delta_a = 0.0481857505705;                         % [1/rad]
aircraft.CY_delta_r = -0.24822;                         % [1/rad]

% Rolling moment coefficient:
aircraft.Cl_beta_0 = -0.164200;                          % [1/rad]
aircraft.Cl_beta_alpha = -3.197;                      % [1/rad^2]
aircraft.Cl_p = -0.595325;                               % [1/rad]
aircraft.Cl_r_0 = 0.101273;                             % [1/rad]
aircraft.Cl_r_alpha = 2.918;                         % [1/rad^2]
aircraft.Cl_delta_a = -0.25965;                         % [1/rad]
aircraft.Cl_delta_r = 0.03180;                         % [1/rad]

% Yawing moment coefficient:
aircraft.Cn_beta = 0.183052;                            % [1/rad]
aircraft.Cn_p_0 = 0.018399;                             % [1/rad]
aircraft.Cn_p_alpha = -1.542;                         % [1/rad^2]
aircraft.Cn_r = -0.228938;                               % [1/rad]
aircraft.Cn_delta_a_0 = -0.00143;                       % [1/rad]
aircraft.Cn_delta_a_alpha = -0.190;                   % [1/rad^2]
aircraft.Cn_delta_r = -0.12800;                         % [1/rad]
