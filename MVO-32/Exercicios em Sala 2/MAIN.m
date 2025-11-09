clear
close all
clc

global g
global aircraft

g = 9.80665;

aircraft = struct('S',260,'c',6.6,'b',44.8,...
    'm',130000,'Ixx',6.011e6,'Iyy',10.53e6,...
    'Izz',15.73e6,'Ixz',0.33e6,...
    'i_p_deg',2.17,'x_p',7.50,'z_p',2.65,...
    'Tmax',452000,'n_rho',0.8);

psidot_deg_s_eq = 0; % com essa variável sendo zero o voo é reto e nivelado

trim_par(1) = struct('V',77,'h',600,'gamma_deg',-3,...
    'thetadot_deg_s',0,'psidot_deg_s',psidot_deg_s_eq);
trim_par(2) = struct('V',131.5,'h',3000,'gamma_deg',0,...
    'thetadot_deg_s',0,'psidot_deg_s',psidot_deg_s_eq);
trim_par(3) = struct('V',264,'h',10000,'gamma_deg',0,...
    'thetadot_deg_s',0,'psidot_deg_s',psidot_deg_s_eq);

% Test case:
Xtest = [
    250
    2
    5
    5
    5000
    0
    2
    10
    10
    15
    10
    0
    ];
Utest = [0.5; -5; 5; 5; 5];
dynamics(0,Xtest,Utest,3)

% https://www.mathworks.com/help/optim/ug/tolerances-and-stopping-criteria.html
options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);
% % In more recent MATLAB versions, one of the following alternatives may be used:
% options = optimoptions(@fsolve,'Display','iter','StepTolerance',1e-10,'FunctionTolerance',1e-10);
% % or:
% options = optimoptions(@fsolve,'Display','iter','TolX',1e-10,'TolFun',1e-10);

%--------------------------------------------------------------------------
% Ex. 2) Trimmed flight conditions:

trim_output(1:3) = struct('X_eq',zeros(12,1),'U_eq',zeros(5,1),...
    'Y_eq',zeros(11,1));

for i_cond=1:3
    x_eq_0 = zeros(12,1);
    x_eq_0(1) = trim_par(i_cond).V;
    x_eq = fsolve(@trim_function,x_eq_0,options,trim_par,i_cond);
    [~,X_eq,U_eq,Y_eq] = trim_function(x_eq,trim_par,i_cond);
    
    trim_output(i_cond).X_eq = X_eq;
    trim_output(i_cond).U_eq = U_eq;
    trim_output(i_cond).Y_eq = Y_eq;

    fprintf('----- A%d FLIGHT CONDITION -----\n\n',i_cond);
    fprintf('   %-10s = %10.4f %-4s\n','gamma',trim_par(i_cond).gamma_deg,'deg');
    fprintf('   %-10s = %10.4f %-4s\n','theta_dot',trim_par(i_cond).thetadot_deg_s,'deg/s');
    fprintf('   %-10s = %10.4f %-4s\n','psi_dot',trim_par(i_cond).psidot_deg_s,'deg/s');
    fprintf('\n');
    fprintf('   %-10s = %10.2f %-4s\n','V',X_eq(1),'m/s');
    fprintf('   %-10s = %10.4f %-4s\n','alpha',X_eq(2),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','q',X_eq(3),'deg/s');
    fprintf('   %-10s = %10.4f %-4s\n','theta',X_eq(4),'deg');
    fprintf('   %-10s = %10.1f %-4s\n','h',X_eq(5),'m');
    fprintf('\n');
    fprintf('   %-10s = %10.4f %-4s\n','beta',X_eq(7),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','phi',X_eq(8),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','p',X_eq(9),'deg/s');
    fprintf('   %-10s = %10.4f %-4s\n','r',X_eq(10),'deg/s');
    fprintf('   %-10s = %10.4f %-4s\n','psi',X_eq(11),'deg');
    fprintf('\n');
    fprintf('   %-10s = %10.2f %-4s\n','throttle',U_eq(1)*100,'%');
    fprintf('   %-10s = %10.2f %-4s\n','Thrust',Y_eq(2),'N');
    fprintf('   %-10s = %10.4f %-4s\n','i_t',U_eq(2),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','delta_e',U_eq(3),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','delta_a',U_eq(4),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','delta_r',U_eq(5),'deg');
    fprintf('\n');
    fprintf('   %-10s = %10.4f %-4s\n','Mach',Y_eq(3),'');
    fprintf('   %-10s = %10.4f %-4s\n','rho',Y_eq(10),'kg/m³');
    fprintf('   %-10s = %10.2f %-4s\n','q bar',Y_eq(11),'N/m²');
    fprintf('\n');
    fprintf('   %-10s = %10.4f %-4s\n','C_D',Y_eq(4),'');
    fprintf('   %-10s = %10.4f %-4s\n','C_L',Y_eq(5),'');
    fprintf('   %-10s = %10.4f %-4s\n','C_m',Y_eq(6),'');
    fprintf('\n');
    fprintf('   %-10s = %10.4f %-4s\n','C_Y',Y_eq(7),'');
    fprintf('   %-10s = %10.4f %-4s\n','C_l',Y_eq(8),'');
    fprintf('   %-10s = %10.4f %-4s\n','C_n',Y_eq(9),'');
    fprintf('\n');
end

save trim_output.mat trim_output


%--------------------------------------------------------------------------
% Ex. 2) Simulation of trimmed flight condition:

% i_cond = 3;
% 
% X0 = trim_output(i_cond).X_eq;
% U0 = trim_output(i_cond).U_eq;
% 
% dt = 0.010;
% t_f = 20;
% T = 0:dt:t_f;
% [X,Y] = ode4xy(@simulate_flight,T,X0,U0,i_cond);
% U = Y(:,end-4:end);
% Y = Y(:,1:end-5);
% 
% plot_long
% plot_latdir
% plot_controls
% plot_path
% plot_outputs


% %--------------------------------------------------------------------------
% % Ex. 3) Simulation of trimmed flight condition:
% 
% i_cond = 3;
% 
% X0 = trim_output(i_cond).X_eq;
% U0 = trim_output(i_cond).U_eq;
% 
% dt = 0.010;
% t_f = 360/psidot_deg_s_eq;
% T = 0:dt:t_f;
% [X,Y] = ode4xy(@simulate_flight,T,X0,U0,i_cond);
% U = Y(:,end-4:end);
% Y = Y(:,1:end-5);
% 
% plot_long
% plot_latdir
% plot_controls
% plot_path
% plot_outputs


% %--------------------------------------------------------------------------
% Ex. 4) Linearization around trimmed flight conditions:

nX = length(X_eq);
nU = length(U_eq);
nY = length(Y_eq);

lin_output(1:3) = struct('A',zeros(nX,nX),...
    'B',zeros(nX,nU),...
    'C',zeros(nY,nX),...
    'D',zeros(nY,nU));

delta_val = 1e-5;

for i_cond=1:3
    X_eq = trim_output(i_cond).X_eq;
    U_eq = trim_output(i_cond).U_eq;

    A = zeros(nX,nX);
    C = zeros(nY,nX);
    for j=1:nX
        dX = zeros(nX,1);
        dX(j) = delta_val;
        [Xdot_plus,Y_plus] = dynamics(0,X_eq+dX,U_eq,i_cond);
        [Xdot_minus,Y_minus] = dynamics(0,X_eq-dX,U_eq,i_cond);
        A(:,j) = (Xdot_plus - Xdot_minus)/(2*dX(j));
        C(:,j) = (Y_plus - Y_minus)/(2*dX(j));
    end

    B = zeros(nX,nU);
    D = zeros(nY,nU);
    for j=1:nU
        dU = zeros(nU,1);
        dU(j) = delta_val;
        [Xdot_plus,Y_plus] = dynamics(0,X_eq,U_eq+dU,i_cond);
        [Xdot_minus,Y_minus] = dynamics(0,X_eq,U_eq-dU,i_cond);
        B(:,j) = (Xdot_plus - Xdot_minus)/(2*dU(j));
        D(:,j) = (Y_plus - Y_minus)/(2*dU(j));
    end

    lin_output(i_cond).A = A;
    lin_output(i_cond).B = B;
    lin_output(i_cond).C = C;
    lin_output(i_cond).D = D;

end

save lin_output.mat lin_output

% [eigvec,eigval] = eig(lin_output(3).A);
sel_states = [1:5 7:11]; % importante fazer isso para tirar x e y da análise dos autovetores, pois eles podem atrapalhar com os valores
[eigvec,eigval] = eig(lin_output(3).A(sel_states,sel_states));
eigval = diag(eigval);
[eigval,i_sort] = sort(eigval);
eigvec = eigvec(:,i_sort);

damp(eigval) % matlab já tem uma forma da calcular damp diretamente


% %--------------------------------------------------------------------------
% % Ex. 5) Simulation of disturbed initial condition:
% 
% i_cond = 3;
% 
% X0 = trim_output(i_cond).X_eq;
% X0(7) = X0(7) + 2;
% U0 = trim_output(i_cond).U_eq;
% 
% dt = 0.010;
% t_f = 20;
% T = 0:dt:t_f;
% [X,Y] = ode4xy(@simulate_flight,T,X0,U0,i_cond);
% U = Y(:,end-4:end);
% Y = Y(:,1:end-5);
% 
% plot_long
% plot_latdir
% plot_controls
% plot_path
% plot_outputs


% %--------------------------------------------------------------------------
% % Ex. 6) Simulation of response to a rudder doublet:
% 
% i_cond = 3;
% 
% X0 = trim_output(i_cond).X_eq;
% U0 = trim_output(i_cond).U_eq;
% 
% dt = 0.010;
% t_f = 30;
% T = 0:dt:t_f;
% [X,Y] = ode4xy(@simulate_flight_Ex6,T,X0,U0,i_cond);
% U = Y(:,end-4:end);
% Y = Y(:,1:end-5);
% 
% plot_long
% plot_latdir
% plot_controls
% plot_path
% plot_outputs


% %--------------------------------------------------------------------------
% Ex. 7) Simulation of disturbed initial condition:

i_cond = 3;

X0 = trim_output(i_cond).X_eq;
X0(8) = X0(8) + 20;
U0 = trim_output(i_cond).U_eq;

dt = 0.010;
t_f = 500;
T = 0:dt:t_f;
[X,Y] = ode4xy(@simulate_flight,T,X0,U0,i_cond);
U = Y(:,end-4:end);
Y = Y(:,1:end-5);

plot_long
plot_latdir
plot_controls
plot_path
plot_outputs