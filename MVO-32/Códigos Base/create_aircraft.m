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
aircraft.Cm_i_t = -0.080744;          % [1/rad] (Cmg2 - AVL ST)
aircraft.Cm_delta_e = -0.063284;        % [1/rad] (Cmd3 - AVL ST)
% Drag coefficient:
aircraft.CD_0 = 0.06378;            % [-] (CD0)
aircraft.CD_alpha = 0.002419;             % [1/rad]
aircraft.CD_alpha2 = 0.031170;             % [1/rad^2]
aircraft.CD_q = 0.410019;           % [1/rad] (CXq - AVL SB)
aircraft.CD_i_t_0 = -0.02452;         % [1/rad] (CDffg2 * 180/pi)
aircraft.CD_i_t_alpha = 0.013070;       % [1/rad^2] (CDit,alpha)
aircraft.CD_delta_e_0 = -0.01821;       % [1/rad] (CDffd3 * 180/pi)
aircraft.CD_delta_e_alpha = 0.009777;     % [1/rad^2] (CDdelta_e,alpha)
aircraft.CD_beta2 = 0.201532;         % [1/rad^2] (CD/beta^2)
aircraft.CD_p2 = -0.84032;           % [1/rad^2] (CDp^2)
aircraft.CD_r2 = -0.13295;           % [1/rad^2] (CDr^2)
aircraft.CD_delta_a2 = 0.32960;                 % [1/rad^2] (CDdelta_a^2)
aircraft.CD_delta_r2 = 0.00000;        % [1/rad^2] (CDdelta_r^2)
% Side force coefficient:
aircraft.CY_beta = -0.527902;         % [1/rad] (CYb - AVL ST)
aircraft.CY_p_0 = 0.143403;          % [1/rad] (CYp * -1 for alpha=0)
aircraft.CY_p_alpha = -1.966691;        % [1/rad^2] (CYp,alpha)
aircraft.CY_r = -0.046680;           % [1/rad] (CYr - AVL ST)
aircraft.CY_delta_a = 0.000000;        % [1/rad] (CYd2 - NOVO)
aircraft.CY_delta_r = 0.18320;         % [1/rad] (CYd4 * -1 * 180/pi)
% Rolling moment coefficient:
aircraft.Cl_beta_0 = -0.090725;        % [1/rad] (Clb for alpha=0)
aircraft.Cl_beta_alpha = -0.901288;      % [1/rad^2] (Clbeta,alpha)
aircraft.Cl_p = -0.600989;           % [1/rad] (Clp - AVL ST)
aircraft.Cl_r_0 = 0.075773;          % [1/rad] (Clr for alpha=0)
aircraft.Cl_r_alpha = -3.993998;        % [1/rad^2] (Clr,alpha)
aircraft.Cl_delta_a = 0.000000;        % [1/rad] (Cld2)
aircraft.Cl_delta_r = 0.02269;         % [1/rad] (Cld4 * -1 * 180/pi)
% Yawing moment coefficient:
aircraft.Cn_beta = -0.068256;         % [1/rad] (Cnb - AVL ST)
aircraft.Cn_p_0 = 0.013035;          % [1/rad] (Cnp for alpha=0)
aircraft.Cn_p_alpha = -0.955782;        % [1/rad^2] (Cnp,alpha)
aircraft.Cn_r = -0.153312;           % [1/rad] (Cnr - AVL ST)
aircraft.Cn_delta_a_0 = 0.000000;       % [1/rad] (Cnd2 * -1 * 180/pi for alpha=0)
aircraft.Cn_delta_a_alpha = 0.000000;     % [1/rad^2] (Cndelta_a,alpha)
aircraft.Cn_delta_r = -0.08404;        % [1/rad] (Cnd4 * -1 * 180/pi)
