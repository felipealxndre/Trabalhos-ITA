figure(3)
subplot(231)
plot(T_Ex5,Y_Ex5(:,1))
hold all
grid on
xlabel('t [s]')
ylabel('\gamma [deg]')

subplot(232)
plot(T_Ex5,Y_Ex5(:,8))
hold all
grid on
xlabel('t [s]')
ylabel('Dyn. press. [N/m²]')

subplot(233)
plot(T_Ex5,Y_Ex5(:,3))
hold all
grid on
xlabel('t [s]')
ylabel('Mach')

subplot(234)
plot(T_Ex5,Y_Ex5(:,4))
hold all
grid on
xlabel('t [s]')
ylabel('C_D')

subplot(235)
plot(T_Ex5,Y_Ex5(:,5))
hold all
grid on
xlabel('t [s]')
ylabel('C_L')

subplot(236)
plot(T_Ex5,Y_Ex5(:,6))
hold all
grid on
xlabel('t [s]')
ylabel('C_m')
