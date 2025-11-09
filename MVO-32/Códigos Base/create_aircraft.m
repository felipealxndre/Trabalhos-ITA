% This script creates the structure variable 'aircraft'

% DO NOT change the fields' names!

%--------------------------------------------------------------------------
% Geometric properties:
% This script creates the structure variable 'aircraft'
% DO NOT change the fields' names!
%--------------------------------------------------------------------------
% Geometric properties:
% Valores estimados com base no arquivo .avl
aircraft.S = 105.5;     % Reference wing area [m^2]
aircraft.c = 3.526898620375959;     % Wing mean aerodynamic chord [m]
aircraft.b = 32.803963175201865;     % Wing span [m]
%--------------------------------------------------------------------------
% Inertia properties:
aircraft.m = 44278.940051969854;       % Mass [kg]
aircraft.Ixx = 512653.72254536394;      % I_{xx} [kg.m^2]
aircraft.Iyy = 2389904.9749673614;      % I_{yy} [kg.m^2]
aircraft.Izz = 2748106.383976222;      % I_{zz} [kg.m^2]
aircraft.Ixz = 39043.319615114786;      % I_{xz} [kg.m^2]
%--------------------------------------------------------------------------
% Propulsive model properties:
aircraft.iota_1_deg = 2;  % Left engine incidence [deg]
aircraft.tau_1_deg = 1.5;   % Left engine toe-in angle [deg]
aircraft.x_1 = 12.0;      % Left engine thrust action point x_b coordinate [m]
aircraft.y_1 = -4.2;      % Left engine thrust action point y_b coordinate [m]
aircraft.z_1 = -2.5192;      % Left engine thrust action point z_b coordinate [m]

aircraft.iota_2_deg = 2;  % Right engine incidence [deg]
aircraft.tau_2_deg = -1.5;   % Right engine toe-in angle [deg]
aircraft.x_2 = 12.0;      % Right engine thrust action point x_b coordinate [m]
aircraft.y_2 = 4.2;      % Right engine thrust action point y_b coordinate [m]
aircraft.z_2 = -2.5192;      % Right engine thrust action point z_b coordinate [m]


aircraft.Tmax = 100000;     % Maximum sea-level thrust for a single engine [N]
aircraft.n_rho = 0.8;     % Air density exponent in the thrust model [-]
%--------------------------------------------------------------------------
% Aerodynamic data bank:
% Lift coefficient:
aircraft.CL_0 = 0.09059;            % [-] (CL0)
aircraft.CL_alpha = 7.357614;         % [1/rad] (CLa - AVL ST)
aircraft.CL_q = 17.315289;           % [1/rad] (CLq - AVL ST)
aircraft.CL_i_t = 0.021783;          % [1/rad] (CLg2 - AVL ST)
aircraft.CL_delta_e = 0.016056;        % [1/rad] (CLd3 - AVL ST)
% Pitching moment coefficient:
aircraft.Cm_0 = -0.10796;           % [-] (Cm0)
aircraft.Cm_alpha = -3.193383;         % [1/rad] (Cma - AVL ST)
aircraft.Cm_q = -48.697758;          % [1/rad] (Cmq - AVL ST)
aircraft.Cm_i_t = -4.62643;          % [1/rad] (Cmg2 - AVL ST)
aircraft.Cm_delta_e = -3.62604;        % [1/rad] (Cmd3 - AVL ST)
% Drag coefficient:
aircraft.CD_0 = 0.06378;            % [-] (CD0)
aircraft.CD_alpha = 0.002419;             % [1/rad]
aircraft.CD_alpha2 = 0.031170;             % [1/rad^2]
aircraft.CD_q = 0.410019;           % [1/rad] (CXq - AVL SB)
aircraft.CD_i_t_0 = 0.01914;         % [1/rad] (CDffg2 * 180/pi)
aircraft.CD_i_t_alpha = 0.013070;       % [1/rad^2] (CDit,alpha)
aircraft.CD_delta_e_0 = 0.01444;       % [1/rad] (CDffd3 * 180/pi)
aircraft.CD_delta_e_alpha = 0.009777;     % [1/rad^2] (CDdelta_e,alpha)
aircraft.CD_beta2 = -0.26262;         % [1/rad^2] (CD/beta^2)
aircraft.CD_p2 = -0.84305;           % [1/rad^2] (CDp^2)
aircraft.CD_r2 = -0.13242;           % [1/rad^2] (CDr^2)
aircraft.CD_delta_a2 = 0.32960;                 % [1/rad^2] (CDdelta_a^2)
aircraft.CD_delta_r2 = 0.00000;        % [1/rad^2] (CDdelta_r^2)
% Side force coefficient:
aircraft.CY_beta = 0.527902;         % [1/rad] (CYb - AVL ST)
aircraft.CY_p_0 = 0.025974;          % [1/rad] (CYp * -1 for alpha=0)
aircraft.CY_p_alpha = -1.966691;        % [1/rad^2] (CYp,alpha)
aircraft.CY_r = 0.048275;           % [1/rad] (CYr - AVL ST)
aircraft.CY_delta_a = 0.000000;        % [1/rad] (CYd2 - NOVO)
aircraft.CY_delta_r = -0.18320;         % [1/rad] (CYd4 * -1 * 180/pi)
% Rolling moment coefficient:
aircraft.Cl_beta_0 = -0.143271;        % [1/rad] (Clb for alpha=0)
aircraft.Cl_beta_alpha = -0.901288;      % [1/rad^2] (Clbeta,alpha)
aircraft.Cl_p = -0.606123;           % [1/rad] (Clp - AVL ST)
aircraft.Cl_r_0 = 0.130652;          % [1/rad] (Clr for alpha=0)
aircraft.Cl_r_alpha = -3.993998;        % [1/rad^2] (Clr,alpha)
aircraft.Cl_delta_a = 0.000000;        % [1/rad] (Cld2)
aircraft.Cl_delta_r = 0.02269;         % [1/rad] (Cld4 * -1 * 180/pi)
% Yawing moment coefficient:
aircraft.Cn_beta = -0.068256;         % [1/rad] (Cnb - AVL ST)
aircraft.Cn_p_0 = -0.069117;          % [1/rad] (Cnp for alpha=0)
aircraft.Cn_p_alpha = -0.955782;        % [1/rad^2] (Cnp,alpha)
aircraft.Cn_r = -0.153312;           % [1/rad] (Cnr - AVL ST)
aircraft.Cn_delta_a_0 = 0.000000;       % [1/rad] (Cnd2 * -1 * 180/pi for alpha=0)
aircraft.Cn_delta_a_alpha = 0.000000;     % [1/rad^2] (Cndelta_a,alpha)
aircraft.Cn_delta_r = -0.08417;        % [1/rad] (Cnd4 * -1 * 180/pi)


% 
% % This script creates the structure variable 'aircraft'
% 
% % DO NOT change the fields' names!
% 
% %--------------------------------------------------------------------------
% % Geometric properties:
% aircraft.S = 81.25;              % Reference wing area [m^2] #ok!
% aircraft.c = 2.53;               % Wing mean aerodynamic chord [m] #ok!
% aircraft.b = 35.99;              % Wing span [m] #ok!
% 
% %--------------------------------------------------------------------------
% % Inertia properties:
% aircraft.m = 41347.03;           % Mass [kg] #ok!
% aircraft.Ixx = 474148.04;        % I_{xx} [kg.m^2] #ok!
% aircraft.Iyy = 2527894.43;       % I_{yy} [kg.m^2] #ok!
% aircraft.Izz = 2948507.61;       % I_{zz} [kg.m^2] #ok!
% aircraft.Ixz = 19981.74;         % I_{xz} [kg.m^2] #ok!
% 
% %--------------------------------------------------------------------------
% % Propulsive model properties:
% aircraft.iota_1_deg = 2.00;      % Left engine incidence [deg] #ok!
% aircraft.tau_1_deg = 1.50;       % Left engine toe-in angle [deg] #ok!
% aircraft.x_1 = 13.90;            % Left engine thrust action point x_b coordinate [m] #ok!
% aircraft.y_1 = -5.00;            % Left engine thrust action point y_b coordinate [m] #ok!
% aircraft.z_1 = 1.00;             % Left engine thrust action point z_b coordinate [m] #ok!
% 
% aircraft.iota_2_deg = 2.00;      % Right engine incidence [deg] #ok!
% aircraft.tau_2_deg = -1.50;       % Right engine toe-in angle [deg] #ok!
% aircraft.x_2 = 13.90;            % Right engine thrust action point x_b coordinate [m] #ok!
% aircraft.y_2 = 5.00;             % Right engine thrust action point y_b coordinate [m] #ok!
% aircraft.z_2 = 1.00;             % Right engine thrust action point z_b coordinate [m] #ok!
% 
% aircraft.Tmax = 116445.95;       % Maximum sea-level thrust for a single engine [N] #ok!
% aircraft.n_rho = 0.8;            % Air density exponent in the thrust model [-] #ok!
% 
% %--------------------------------------------------------------------------
% % Aerodynamic data bank:
% 
% % Important notes:
% % -> Aerodynamic moments calculated about the nominal CG, in the
% % directions of the body axes
% % -> Angular velocity components refer to rotational motion about the
% % nominal CG
% % -> All angles used in the aerodynamic databank will need to be converted 
% % to radians
% % -> All angular velocities used in the aerodynamic databank will need to 
% % be converted to radians per second and then made dimensionless by the
% % appropriate multiplication by c/(2*V) or b/(2*V)
% 
% % Lift coefficient:
% aircraft.CL_0 = -0.0407;                             % [-] #ok!
% aircraft.CL_alpha = 8.44;                            % [1/rad] #ok!
% aircraft.CL_q = 19.4;                                % [1/rad] #ok!
% aircraft.CL_i_t = 1.33;                              % [1/rad] #ok!
% aircraft.CL_delta_e = 0.911;                         % [1/rad] #ok!
% 
% % Pitching moment coefficient:
% aircraft.Cm_0 = -0.0596;                             % [-] #ok!
% aircraft.Cm_alpha = -3.06;                           % [1/rad] #ok!
% aircraft.Cm_q = -110.0;                              % [1/rad] #ok!
% aircraft.Cm_i_t = -7.82;                             % [1/rad] #ok!
% aircraft.Cm_delta_e = -5.60;                         % [1/rad] #ok!
% 
% % Drag coefficient:
% aircraft.CD_0 = 0.0190;                              % [-] #ok!
% aircraft.CD_alpha = 0.0;                             % [1/rad] #ok!
% aircraft.CD_alpha2 = 0.0005;                         % [1/rad^2] #ok!
% aircraft.CD_q = 0.516;                               % [1/rad] #ok!
% aircraft.CD_i_t_0 = -0.0256;                         % [1/rad] #ok!
% aircraft.CD_i_t_alpha = 0.3914;                      % [1/rad^2] #ok!
% aircraft.CD_delta_e_0 = -0.0175;                     % [1/rad] #ok!
% aircraft.CD_delta_e_alpha = 0.2677;                  % [1/rad^2] #ok!
% aircraft.CD_beta2 = -0.55282;                         % [1/rad^2] #ok!
% aircraft.CD_p2 = -1.244;                             % [1/rad^2] #ok!
% aircraft.CD_r2 = -0.432;                             % [1/rad^2] #ok!
% aircraft.CD_delta_a2 = 0.2284;                       % [1/rad^2] #ok!
% aircraft.CD_delta_r2 = 0.1234;                       % [1/rad^2] #ok!
% 
% % Side force coefficient:
% aircraft.CY_beta = 0.857;                            % [1/rad] #ok!
% aircraft.CY_p_0 = 0.226;                             % [1/rad] #ok!
% aircraft.CY_p_alpha = -4.9023;                       % [1/rad^2] #ok!
% aircraft.CY_r = -0.0799;                             % [1/rad] #ok!
% aircraft.CY_delta_a = 0.1197;                           % [1/rad] #ok! (AVL)
% aircraft.CY_delta_r = -0.192;                        % [1/rad] #ok!
% 
% % Rolling moment coefficient:
% aircraft.Cl_beta_0 = -2.00e-6;                       % [1/rad] #ok!
% aircraft.Cl_beta_alpha = 3.0e-5;                     % [1/rad^2] #ok!
% aircraft.Cl_p = -0.716;                              % [1/rad] #ok!
% aircraft.Cl_r_0 = 0.0618;                            % [1/rad] #ok!
% aircraft.Cl_r_alpha = 0.9504;                        % [1/rad^2] #ok!
% aircraft.Cl_delta_a = -0.550;                        % [1/rad] #ok!
% aircraft.Cl_delta_r = 0.0215;                        % [1/rad] #ok!
% 
% % Yawing moment coefficient:
% aircraft.Cn_beta = 0.00445;                         % [1/rad] #ok! Sinal contrario do avl
% aircraft.Cn_p_0 = 0.0289;                            % [1/rad] #ok!
% aircraft.Cn_p_alpha = -1.5474;                       % [1/rad^2] #ok!
% aircraft.Cn_r = -0.223;                              % [1/rad] #ok!
% aircraft.Cn_delta_a_0 = -0.012;                      % [1/rad] #ok!
% aircraft.Cn_delta_a_alpha = 0.1881;                  % [1/rad^2] #ok!
% aircraft.Cn_delta_r = -0.0939;                       % [1/rad] #ok!
% 