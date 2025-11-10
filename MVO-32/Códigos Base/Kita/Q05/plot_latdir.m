figure(2)

subplot(231)
plot(T,X(:,7))
hold all
grid on
xlabel('t [s]')
ylabel('\beta [deg]')

subplot(232)
plot(T,X(:,8))
hold all
grid on
xlabel('t [s]')
ylabel('\phi [deg]')

subplot(233)
plot(T,X(:,9))
hold all
grid on
xlabel('t [s]')
ylabel('p [deg/s]')

subplot(234)
plot(T,X(:,10))
hold all
grid on
xlabel('t [s]')
ylabel('r [deg/s]')

subplot(235)
plot(T,X(:,11))
hold all
grid on
xlabel('t [s]')
ylabel('\psi [deg]')

subplot(236)
plot(T,X(:,12))
hold all
grid on
xlabel('t [s]')
ylabel('y [m]')
