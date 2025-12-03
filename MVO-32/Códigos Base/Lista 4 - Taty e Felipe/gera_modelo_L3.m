function [A_L3,B_L3,C_L3,D_L3,X_eq_L3,U_eq_L3,Y_eq_L3] = gera_modelo_L3()
% Gera o ponto de operação e o modelo linearizado da Lista 3 (Ex.2 e Ex.3)
% para a condição de cruzeiro usada na L4.

    clearvars -except g aircraft
    clc

    global g
    global aircraft
    create_aircraft;

    g = 9.80665;

    % Condição de voo A3 (igual à que você usou na MAIN)
    trim_par = struct('V',228.8138886191133,'h',11000);

    options = optimset('Display','off','TolX',1e-10,'TolFun',1e-10);

    % Ex.2: trim de voo reto nivelado em cruzeiro
    x_eq_0 = zeros(5,1);
    x_eq   = fsolve(@trim_function_Ex2,x_eq_0,options,trim_par);
    [~,X_eq_L3,U_eq_L3,Y_eq_L3] = trim_function_Ex2(x_eq,trim_par);

    % Ex.3: linearização numérica em torno do equilíbrio
    nX = length(X_eq_L3);
    nU = length(U_eq_L3);
    nY = length(Y_eq_L3);

    delta_val = 1e-5;

    A_L3 = zeros(nX,nX);
    C_L3 = zeros(nY,nX);
    for j = 1:nX
        dX = zeros(nX,1);
        dX(j) = delta_val;
        [Xdot_plus, Y_plus]  = dynamics(0, X_eq_L3 + dX, U_eq_L3);
        [Xdot_minus, Y_minus]= dynamics(0, X_eq_L3 - dX, U_eq_L3);
        A_L3(:, j) = (Xdot_plus - Xdot_minus)/(2*delta_val);
        C_L3(:, j) = (Y_plus - Y_minus)/(2*delta_val);
    end

    B_L3 = zeros(nX,nU);
    D_L3 = zeros(nY,nU);
    for j = 1:nU
        dU = zeros(nU,1);
        dU(j) = delta_val;
        [Xdot_plus, Y_plus]  = dynamics(0, X_eq_L3, U_eq_L3 + dU);
        [Xdot_minus, Y_minus]= dynamics(0, X_eq_L3, U_eq_L3 - dU);
        B_L3(:, j) = (Xdot_plus - Xdot_minus)/(2*delta_val);
        D_L3(:, j) = (Y_plus - Y_minus)/(2*delta_val);
    end
end
