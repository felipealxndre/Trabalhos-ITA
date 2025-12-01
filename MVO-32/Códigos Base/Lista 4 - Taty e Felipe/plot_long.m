figure(1)
subplot(231)
plot(T_Ex5,X_Ex5(:,1))
hold all
grid on
xlabel('t [s]')
ylabel('V [m/s]')

subplot(232)
plot(T_Ex5,X_Ex5(:,2))
hold all
grid on
xlabel('t [s]')
ylabel('\alpha [deg]')

subplot(233)
plot(T_Ex5,X_Ex5(:,3))
hold all
grid on
xlabel('t [s]')
ylabel('q [deg/s]')

subplot(234)
plot(T_Ex5,X_Ex5(:,4))
hold all
grid on
xlabel('t [s]')
ylabel('\theta [deg]')

subplot(235)
plot(T_Ex5,X_Ex5(:,5))
hold all
grid on
xlabel('t [s]')
ylabel('h [m]')

subplot(236)
plot(T_Ex5,X_Ex5(:,6))
hold all
grid on
xlabel('t [s]')
ylabel('x [m]')
