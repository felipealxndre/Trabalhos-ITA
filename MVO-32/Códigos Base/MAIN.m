clear
close all
clc

global g
global aircraft
create_aircraft;

g = 9.80665;

trim_par = struct('V',228.8138886191133,'h',11000,'gamma_deg',0,...
    'thetadot_deg_s',0);

% https://www.mathworks.com/help/optim/ug/tolerances-and-stopping-criteria.html
options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);
% % In newer MATLAB versions, the following command might be necessary:
% options = optimoptions(@fsolve,'Display','iter','StepTolerance',1e-10,'FunctionTolerance',1e-10);
% % or:
% options = optimoptions(@fsolve,'Display','iter','TolX',1e-10,'TolFun',1e-10);

%% Ex.2

x_eq_0 = zeros(6,1);   
x_eq = fsolve(@trim_function_Ex2,x_eq_0,options,trim_par);
[~,X_eq,U_eq,Y_eq] = trim_function_Ex2(x_eq,trim_par);

fprintf('\n----- EX.2 -----\n\n');
fprintf('   %-12s = %8.2f %s\n','V',X_eq(1),'m/s');
fprintf('   %-12s = %8.3f %s\n','alpha',X_eq(2),'deg');
fprintf('   %-12s = %8.3f %s\n','q',X_eq(3),'deg/s');
fprintf('   %-12s = %8.3f %s\n','theta',X_eq(4),'deg');
fprintf('   %-12s = %8.1f %s\n','h',X_eq(5),'m');
fprintf('   %-12s = %8.3f %s\n','beta',X_eq(7),'deg');
fprintf('   %-12s = %8.3f %s\n','phi',X_eq(8),'deg');
fprintf('   %-12s = %8.3f %s\n','p',X_eq(9),'deg/s');
fprintf('   %-12s = %8.3f %s\n','r',X_eq(10),'deg/s');
fprintf('\n');
fprintf('   %-12s = %8.2f %s\n','delta1',U_eq(1),'%');
fprintf('   %-12s = %8.2f %s\n','delta2',U_eq(2),'%');
fprintf('   %-12s = %8.3f %s\n','i_t',U_eq(3),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_e',U_eq(4),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_a',U_eq(5),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_r',U_eq(6),'deg');
fprintf('\n');
fprintf('   %-12s = %8.3f %s\n','gamma',Y_eq(1),'deg');
fprintf('   %-12s = %8.1f %s\n','T1',Y_eq(2),'N');
fprintf('   %-12s = %8.1f %s\n','T2',Y_eq(3),'N');
fprintf('   %-12s = %8.3f %s\n','Mach',Y_eq(4),'');
fprintf('   %-12s = %8.4f %s\n','C_D',Y_eq(5),'');
fprintf('   %-12s = %8.4f %s\n','C_L',Y_eq(6),'');
fprintf('   %-12s = %8.4f %s\n','C_m',Y_eq(7),'');
fprintf('   %-12s = %8.4f %s\n','C_Y',Y_eq(8),'');
fprintf('   %-12s = %8.4f %s\n','C_l',Y_eq(9),'');
fprintf('   %-12s = %8.4f %s\n','C_n',Y_eq(10),'');
fprintf('   %-12s = %8.4f %s\n','rho',Y_eq(11),'kg/m^3');
fprintf('   %-12s = %8.1f %s\n','qbar',Y_eq(12),'N/m^2');
fprintf('\n');

%% Ex.3

nX = length(X_eq);
nU = length(U_eq);
nY = length(Y_eq);

lin_output = struct('A',zeros(nX,nX),...
    'B',zeros(nX,nU),...
    'C',zeros(nY,nX),...
    'D',zeros(nY,nU));

delta_val = 1e-5;

A = zeros(nX,nX);
C = zeros(nY,nX);
for j=1:nX
    dX = zeros(nX,1);
    dX(j) = delta_val;
    [Xdot_plus, Y_plus]= dynamics(0, X_eq + dX, U_eq);
    [Xdot_minus, Y_minus]= dynamics(0, X_eq - dX, U_eq);
    A(:, j) = (Xdot_plus - Xdot_minus)/(2*dX(j));
    C(:, j) = (Y_plus - Y_minus)/(2*dX(j));
end

B = zeros(nX,nU);
D = zeros(nY,nU);
for j=1:nU
    dU = zeros(nU,1);
    dU(j) = delta_val;
    [Xdot_plus, Y_plus]= dynamics(0, X_eq, U_eq + dU);
    [Xdot_minus, Y_minus]= dynamics(0, X_eq, U_eq - dU);
    B(:, j) = (Xdot_plus - Xdot_minus)/(2*dU(j));
    D(:, j) = (Y_plus - Y_minus)/(2*dU(j));
end

lin_output.A = A;
lin_output.B = B;
lin_output.C = C;
lin_output.D = D;

% [eigvec,eigval] = eig(lin_output.A);
sel_states = 1:10;
[eigvec,eigval] = eig(lin_output.A(sel_states,sel_states));
eigval = diag(eigval);
[eigval,i_sort] = sort(eigval);
eigvec = eigvec(:,i_sort);

fprintf('\n----- EX.3 -----\n\n');
damp(eigval)

%% Ex.4

x0 = [ X_eq(2); U_eq(1); U_eq(3); 0.0; 0.5; -1.0 ]; % chute: [alpha, δ1, i_t, δe, δa, δr]
[x_sol, ~] = fsolve(@trim_function_Ex4, x0, options, trim_par);
[~, X_eq_fail, U_eq_fail, Y_eq_fail] = trim_function_Ex4(x_sol, trim_par);

fprintf('\n----- EX.4 -----\n\n');
fprintf('   %-12s = %8.2f %s\n','V',X_eq_fail(1),'m/s');
fprintf('   %-12s = %8.3f %s\n','alpha',X_eq_fail(2),'deg');
fprintf('   %-12s = %8.3f %s\n','q',X_eq_fail(3),'deg/s');
fprintf('   %-12s = %8.3f %s\n','theta',X_eq_fail(4),'deg');
fprintf('   %-12s = %8.1f %s\n','h',X_eq_fail(5),'m');
fprintf('   %-12s = %8.3f %s\n','beta',X_eq_fail(7),'deg');
fprintf('   %-12s = %8.3f %s\n','phi',X_eq_fail(8),'deg');
fprintf('   %-12s = %8.3f %s\n','p',X_eq_fail(9),'deg/s');
fprintf('   %-12s = %8.3f %s\n','r',X_eq_fail(10),'deg/s');
fprintf('\n');
fprintf('   %-12s = %8.2f %s\n','delta1',U_eq_fail(1),'%');
fprintf('   %-12s = %8.2f %s\n','delta2',U_eq_fail(2),'%');
fprintf('   %-12s = %8.3f %s\n','i_t',U_eq_fail(3),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_e',U_eq_fail(4),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_a',U_eq_fail(5),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_r',U_eq_fail(6),'deg');
fprintf('\n');
fprintf('   %-12s = %8.3f %s\n','gamma',Y_eq_fail(1),'deg');
fprintf('   %-12s = %8.1f %s\n','T1',Y_eq_fail(2),'N');
fprintf('   %-12s = %8.1f %s\n','T2',Y_eq_fail(3),'N');
fprintf('   %-12s = %8.3f %s\n','Mach',Y_eq_fail(4),'');
fprintf('   %-12s = %8.4f %s\n','C_D',Y_eq_fail(5),'');
fprintf('   %-12s = %8.4f %s\n','C_L',Y_eq_fail(6),'');
fprintf('   %-12s = %8.4f %s\n','C_m',Y_eq_fail(7),'');
fprintf('   %-12s = %8.4f %s\n','C_Y',Y_eq_fail(8),'');
fprintf('   %-12s = %8.4f %s\n','C_l',Y_eq_fail(9),'');
fprintf('   %-12s = %8.4f %s\n','C_n',Y_eq_fail(10),'');
fprintf('   %-12s = %8.4f %s\n','rho',Y_eq_fail(11),'kg/m^3');
fprintf('   %-12s = %8.1f %s\n','qbar',Y_eq_fail(12),'N/m^2');
fprintf('\n');

%% Ex.5

X0 = X_eq;
U0 = U_eq;

dt = 0.010;
t_f = 30;
T = 0:dt:t_f;
[X,Y] = ode4xy(@simulate_flight,T,X0,U0);
U = Y(:,end-4:end);
Y = Y(:,1:end-5);

plot_long
plot_latdir
plot_controls
plot_path
plot_outputs

%% Ex.6

x_eq_0 = zeros(6,1);   
x_eq = fsolve(@trim_function_Ex6,x_eq_0,options,trim_par);
[~,X_eq,U_eq,Y_eq] = trim_function_Ex6(x_eq,trim_par);

nX = length(X_eq);
nU = length(U_eq);
nY = length(Y_eq);

lin_output = struct('A',zeros(nX,nX),...
    'B',zeros(nX,nU),...
    'C',zeros(nY,nX),...
    'D',zeros(nY,nU));

delta_val = 1e-5;

A = zeros(nX,nX);
C = zeros(nY,nX);
for j=1:nX
    dX = zeros(nX,1);
    dX(j) = delta_val;
    [Xdot_plus, Y_plus]= dynamics_Ex6(0, X_eq + dX, U_eq);
    [Xdot_minus, Y_minus]= dynamics_Ex6(0, X_eq - dX, U_eq);
    A(:, j) = (Xdot_plus - Xdot_minus)/(2*dX(j));
    C(:, j) = (Y_plus - Y_minus)/(2*dX(j));
end

B = zeros(nX,nU);
D = zeros(nY,nU);
for j=1:nU
    dU = zeros(nU,1);
    dU(j) = delta_val;
    [Xdot_plus, Y_plus]= dynamics_Ex6(0, X_eq, U_eq + dU);
    [Xdot_minus, Y_minus]= dynamics_Ex6(0, X_eq, U_eq - dU);
    B(:, j) = (Xdot_plus - Xdot_minus)/(2*dU(j));
    D(:, j) = (Y_plus - Y_minus)/(2*dU(j));
end

lin_output.A = A;
lin_output.B = B;
lin_output.C = C;
lin_output.D = D;

% [eigvec,eigval] = eig(lin_output.A);
sel_states = 1:10;
[eigvec,eigval] = eig(lin_output.A(sel_states,sel_states));
eigval = diag(eigval);
[eigval,i_sort] = sort(eigval);
eigvec = eigvec(:,i_sort);

fprintf('\n----- EX.6 -----\n\n');
damp(eigval)

%% Ex.7

