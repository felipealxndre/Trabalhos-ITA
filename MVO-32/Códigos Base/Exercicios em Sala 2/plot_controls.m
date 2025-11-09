figure(3)

subplot(231)
plot(T,U(:,1)*100)
hold all
grid on
xlabel('t [s]')
ylabel('Throttle [%]')

subplot(232)
plot(T,U(:,2))
hold all
grid on
xlabel('t [s]')
ylabel('i_t [deg]')

subplot(233)
plot(T,U(:,4))
hold all
grid on
xlabel('t [s]')
ylabel('\delta_a [deg]')

subplot(234)
plot(T,Y(:,2))
hold all
grid on
xlabel('t [s]')
ylabel('Thrust [N]')

subplot(235)
plot(T,U(:,3))
hold all
grid on
xlabel('t [s]')
ylabel('\delta_e [deg]')

subplot(236)
plot(T,U(:,5))
hold all
grid on
xlabel('t [s]')
ylabel('\delta_r [deg]')