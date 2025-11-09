clear
close all
clc

global g
global aircraft

g = 9.80665;

aircraft = struct('S',260,'c',6.6,'b',44.8,...
    'm',130000,'Iyy',10.53e6,...
    'i_p_deg',2.17,'x_p',7.50,'z_p',2.65,...
    'Tmax',452000,'n_rho',0.8);

trim_par(1) = struct('V',77,'h',600,'gamma_deg',-3,...
    'thetadot_deg_s',0);
trim_par(2) = struct('V',131.5,'h',3000,'gamma_deg',0,...
    'thetadot_deg_s',0);
trim_par(3) = struct('V',264,'h',10000,'gamma_deg',0,...
    'thetadot_deg_s',0);

X_test = [
    trim_par(3).V
    2
    10
    5
    trim_par(3).h
    0
    ];

U_test = [
    0.5
    -5
    -5
    ];

Xdot_test = long_dynamics(0,X_test,U_test,3);

% https://www.mathworks.com/help/optim/ug/tolerances-and-stopping-criteria.html
options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);
% % In newer MATLAB versions, the following command might be necessary:
% options = optimoptions(@fsolve,'Display','iter','StepTolerance',1e-10,'FunctionTolerance',1e-10);
% % or:
% options = optimoptions(@fsolve,'Display','iter','TolX',1e-10,'TolFun',1e-10);

%--------------------------------------------------------------------------
% Ex. 2) Trimmed flight conditions:

trim_output(1:3) = struct('X_eq',zeros(6,1),'U_eq',zeros(3,1),...
    'Y_eq',zeros(8,1));

for i_cond=1:3
    x_eq_0 = zeros(6,1);
    x_eq_0(1) = trim_par(i_cond).V;
    x_eq = fsolve(@trim_function,x_eq_0,options,trim_par,i_cond);
    [~,X_eq,U_eq,Y_eq] = trim_function(x_eq,trim_par,i_cond);

    trim_output(i_cond).X_eq = X_eq;
    trim_output(i_cond).U_eq = U_eq;
    trim_output(i_cond).Y_eq = Y_eq;

    fprintf('\n----- A%d FLIGHT CONDITION -----\n\n',i_cond);
    fprintf('   %-10s = %10.4f %-4s\n','gamma',trim_par(i_cond).gamma_deg,'deg');
    fprintf('   %-10s = %10.4f %-4s\n','theta_dot',trim_par(i_cond).thetadot_deg_s,'deg/s');
    fprintf('\n');
    fprintf('   %-10s = %10.2f %-4s\n','V',X_eq(1),'m/s');
    fprintf('   %-10s = %10.4f %-4s\n','alpha',X_eq(2),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','q',X_eq(3),'deg/s');
    fprintf('   %-10s = %10.4f %-4s\n','theta',X_eq(4),'deg');
    fprintf('   %-10s = %10.1f %-4s\n','h',X_eq(5),'m');
    fprintf('\n');
    fprintf('   %-10s = %10.2f %-4s\n','throttle',U_eq(1)*100,'%');
    fprintf('   %-10s = %10.2f %-4s\n','Thrust',Y_eq(2),'N');
    fprintf('   %-10s = %10.4f %-4s\n','i_t',U_eq(2),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','delta_e',U_eq(3),'deg');
    fprintf('\n');
    fprintf('   %-10s = %10.4f %-4s\n','C_D',Y_eq(4),'');
    fprintf('   %-10s = %10.4f %-4s\n','C_L',Y_eq(5),'');
    fprintf('   %-10s = %10.4f %-4s\n','C_m',Y_eq(6),'');
    fprintf('\n');
    fprintf('   %-10s = %10.4f %-4s\n','Mach',Y_eq(3),'');
    fprintf('   %-10s = %10.4f %-4s\n','rho',Y_eq(7),'kg/m³');
    fprintf('   %-10s = %10.2f %-4s\n','q bar',Y_eq(8),'N/m²');
    fprintf('\n');
end

save trim_output.mat trim_output

% %--------------------------------------------------------------------------
% % Ex. 2) Simulation of trimmed flight condition:
% 
% i_cond = 1;
% 
% X0 = trim_output(i_cond).X_eq;
% U0 = trim_output(i_cond).U_eq;
% 
% dt = 0.010;
% t_f = 20;
% T = 0:dt:t_f;
% [X,Y] = ode4xy(@simulate_flight,T,X0,U0,i_cond);
% U = Y(:,end-2:end);
% Y = Y(:,1:end-3);
% 
% plot_long
% plot_controls
% plot_path
% plot_outputs


%--------------------------------------------------------------------------
% Ex. 3) Linearization around trimmed flight conditions:

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
        [Xdot_plus, Y_plus]= long_dynamics(0, X_eq + dX, U_eq, i_cond);
        [Xdot_minus, Y_minus]= long_dynamics(0, X_eq - dX, U_eq, i_cond);
        A(:, j) = (Xdot_plus - Xdot_minus)/(2*dX(j));
        C(:, j) = (Y_plus - Y_minus)/(2*dX(j));
    end
    
    B = zeros(nX,nU);
    D = zeros(nY,nU);
    for j=1:nU
        dU = zeros(nU,1);
        dU(j) = delta_val;
        [Xdot_plus, Y_plus]= long_dynamics(0, X_eq, U_eq + dU, i_cond);
        [Xdot_minus, Y_minus]= long_dynamics(0, X_eq, U_eq - dU, i_cond);
        B(:, j) = (Xdot_plus - Xdot_minus)/(2*dU(j));
        D(:, j) = (Y_plus - Y_minus)/(2*dU(j));
    end
    
    lin_output(i_cond).A = A;
    lin_output(i_cond).B = B;
    lin_output(i_cond).C = C;
    lin_output(i_cond).D = D;

end

save lin_output.mat lin_output

% [eigvec,eigval] = eig(lin_output(3).A);
sel_states = 1:5;
[eigvec,eigval] = eig(lin_output(3).A(sel_states,sel_states));
eigval = diag(eigval);
[eigval,i_sort] = sort(eigval);
eigvec = eigvec(:,i_sort);

damp(eigval)


% %--------------------------------------------------------------------------
% % Ex. 4) Simulation of disturbed initial condition:
% 
% i_cond = 3;
% 
% X0 = trim_output(i_cond).X_eq;
% X0(2) = X0(2) + 2;
% U0 = trim_output(i_cond).U_eq;
% 
% dt = 0.010;
% %t_f = 10;
% t_f = 500;
% T = 0:dt:t_f;
% [X,Y] = ode4xy(@simulate_flight,T,X0,U0,i_cond);
% U = Y(:,end-2:end);
% Y = Y(:,1:end-3);
% 
% plot_long
% plot_controls
% plot_path
% plot_outputs


% %--------------------------------------------------------------------------
% Ex. 5) Simulation of response to an elevator doublet:

%Sem atuador
i_cond = 3;

X0 = trim_output(i_cond).X_eq;
U0 = trim_output(i_cond).U_eq;

dt = 0.010;
t_f = 10;
% t_f = 500;
T = 0:dt:t_f;
[X,Y] = ode4xy(@simulate_flight_Ex5,T,X0,U0,i_cond);
U = Y(:,end-2:end);
Y = Y(:,1:end-3);

plot_long
plot_controls
plot_path
plot_outputs


%Com atuador
i_cond = 3;

X0 = [
     trim_output(i_cond).X_eq
     trim_output(i_cond).U_eq(3)
     ];
U0 = trim_output(i_cond).U_eq;

dt = 0.010;
t_f = 10;
% t_f = 500;
T = 0:dt:t_f;
[X,Y] = ode4xy(@simulate_flight_Ex5_actuator,T,X0,U0,i_cond);
u_e_deg = Y(:, end);
U = Y(:,end-3:end-1);
Y = Y(:,1:end-4);

plot_long
legend('No actuator', 'With actuator')
plot_controls
legend('No actuator', 'With actuator')
plot_path
legend('No actuator', 'With actuator')
plot_outputs
legend('No actuator', 'With actuator')

figure
plot(T, u_e_deg)
hold all
plot(T, X(:,7))
grid on
xlabel('t [s]')
ylabel('Elevator [deg]')
legend('u_e', '\delta_e')

%Dinamica linearizada com atuador
tau_e = 0.050;

for i_cond=1:3
    X_eq = trim_output(i_cond).X_eq;
    U_eq = trim_output(i_cond).U_eq;
    
    A = lin_output(i_cond).A;
    B = lin_output(i_cond).B;
    C = lin_output(i_cond).C;
    D = lin_output(i_cond).D;
    
    A_a = [
        A B(:,3)
        zeros(1,size(A,2)) -1/tau_e
        ];

    B_a =[
        B(:,1:2) zeros(size(B,1), 1)
        zeros(1,2) 1/tau_e
        ];
    C_a = [
        C D(:,3)    
        ];
    D_a = [
        D(:, 1:2) zeros(size(D,1), 1)
        ];

    lin_output(i_cond).A_a = A_a;
    lin_output(i_cond).B_a = B_a;
    lin_output(i_cond).C_a = C_a;
    lin_output(i_cond).D_a = D_a;
end

sel_states = [1:5 7];
[eigvec_a,eigval_a] = eig(lin_output(3).A_a(sel_states,sel_states));
eigval_a = diag(eigval_a);
[eigval_a,i_sort] = sort(eigval_a);
eigvec_a = eigvec_a(:,i_sort);

damp(eigval_a)