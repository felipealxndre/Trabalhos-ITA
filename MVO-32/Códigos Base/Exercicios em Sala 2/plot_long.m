figure(1)

subplot(231)
plot(T,X(:,1))
hold all
grid on
xlabel('t [s]')
ylabel('V [m/s]')

subplot(232)
plot(T,X(:,2))
hold all
grid on
xlabel('t [s]')
ylabel('\alpha [deg]')

subplot(233)
plot(T,X(:,3))
hold all
grid on
xlabel('t [s]')
ylabel('q [deg/s]')

subplot(234)
plot(T,X(:,4))
hold all
grid on
xlabel('t [s]')
ylabel('\theta [deg]')

subplot(235)
plot(T,X(:,5))
hold all
grid on
xlabel('t [s]')
ylabel('h [m]')

subplot(236)
plot(T,X(:,6))
hold all
grid on
xlabel('t [s]')
ylabel('x [m]')