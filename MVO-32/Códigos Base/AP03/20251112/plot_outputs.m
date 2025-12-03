figure(5)

subplot(331)
plot(T,Y(:,1))
hold all
grid on
xlabel('t [s]')
ylabel('\gamma [deg]')

subplot(332)
plot(T,Y(:,11))
hold all
grid on
xlabel('t [s]')
ylabel('Dyn. press. [N/m²]')

subplot(333)
plot(T,Y(:,3))
hold all
grid on
xlabel('t [s]')
ylabel('Mach')

subplot(334)
plot(T,Y(:,4))
hold all
grid on
xlabel('t [s]')
ylabel('C_D')

subplot(335)
plot(T,Y(:,5))
hold all
grid on
xlabel('t [s]')
ylabel('C_L')

subplot(336)
plot(T,Y(:,6))
hold all
grid on
xlabel('t [s]')
ylabel('C_m')

subplot(337)
plot(T,Y(:,7))
hold all
grid on
xlabel('t [s]')
ylabel('C_Y')

subplot(338)
plot(T,Y(:,8))
hold all
grid on
xlabel('t [s]')
ylabel('C_l')

subplot(339)
plot(T,Y(:,9))
hold all
grid on
xlabel('t [s]')
ylabel('C_n')
