figure(2)
subplot(221)
plot(T_Ex5,U_Ex5(:,1)*100)
hold all
grid on
xlabel('t [s]')
ylabel('throttle [%]')

subplot(222)
plot(T_Ex5,U_Ex5(:,2))
hold all
grid on
xlabel('t [s]')
ylabel('i_t [deg]')

subplot(223)
plot(T_Ex5,Y_Ex5(:,2))
hold all
grid on
xlabel('t [s]')
ylabel('Thrust [N]')

subplot(224)
plot(T_Ex5,U_Ex5(:,3))
hold all
grid on
xlabel('t [s]')
ylabel('\delta_e [deg]')
