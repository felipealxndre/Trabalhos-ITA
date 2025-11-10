figure(2)

subplot(231)
plot(T_Ex5,X_Ex5(:,7))
hold all
grid on
xlabel('t [s]')
ylabel('\beta [deg]')

subplot(232)
plot(T_Ex5,X_Ex5(:,8))
hold all
grid on
xlabel('t [s]')
ylabel('\phi [deg]')

subplot(233)
plot(T_Ex5,X_Ex5(:,9))
hold all
grid on
xlabel('t [s]')
ylabel('p [deg/s]')

subplot(234)
plot(T_Ex5,X_Ex5(:,10))
hold all
grid on
xlabel('t [s]')
ylabel('r [deg/s]')

subplot(235)
plot(T_Ex5,X_Ex5(:,11))
hold all
grid on
xlabel('t [s]')
ylabel('\psi [deg]')

subplot(236)
plot(T_Ex5,X_Ex5(:,12))
hold all
grid on
xlabel('t [s]')
ylabel('y [m]')