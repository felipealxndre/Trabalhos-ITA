
x0 = [2;5;0;0;0;0;0];
dt = 0.010;
tF = 20;

% 1D feedback:
Ksim = K_1D;
Ac = A - B*Ksim*C_1D;
Csim = [C_r_w; -Ksim*C_1D];
Dsim = [zeros(1,2); eye(2)];
sys_SAS_1D_sim = ss(Ac,B,Csim,Dsim);
[ysas_1D,tsas_1D,xsas_1D] = initial(sys_SAS_1D_sim,x0,0:dt:tF);

% 1E feedback:
Ksim = K_1E;
Ac = A - B*Ksim*C_1E;
Csim = [C_r_w; -Ksim*C_1E];
Dsim = [zeros(1,2); eye(2)];
sys_SAS_1E_sim = ss(Ac,B,Csim,Dsim);
[ysas_1E,tsas_1E,xsas_1E] = initial(sys_SAS_1E_sim,x0,0:dt:tF);

% No feedback:
Csim = [C_r_w; zeros(2,size(A,1))];
Dsim = [zeros(1,2); eye(2)];
sys_SASoff_sim = ss(A,B,Csim,Dsim);
[y,t,x] = initial(sys_SASoff_sim,x0,0:dt:tF);


figure
subplot(321)
plot(t,x(:,1),'LineWidth',2)
hold all
plot(tsas_1D,xsas_1D(:,1),'LineWidth',2)
plot(tsas_1E,xsas_1E(:,1),'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('\Delta\beta [deg]')

subplot(322)
plot(t,x(:,2),'LineWidth',2)
hold all
plot(tsas_1D,xsas_1D(:,2),'LineWidth',2)
plot(tsas_1E,xsas_1E(:,2),'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('\Delta\phi [deg]')

subplot(323)
plot(t,x(:,3),'LineWidth',2)
hold all
plot(tsas_1D,xsas_1D(:,3),'LineWidth',2)
plot(tsas_1E,xsas_1E(:,3),'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('p [deg/s]')

subplot(324)
plot(t,x(:,4),'LineWidth',2)
hold all
plot(tsas_1D,xsas_1D(:,4),'LineWidth',2)
plot(tsas_1E,xsas_1E(:,4),'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('r [deg/s]')

subplot(325)
plot(t,x(:,5),'LineWidth',2)
hold all
plot(tsas_1D,xsas_1D(:,5),'LineWidth',2)
plot(tsas_1E,xsas_1E(:,5),'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('\Delta{\delta_a} [deg]')

subplot(326)
plot(t,x(:,6),'LineWidth',2)
hold all
plot(tsas_1D,xsas_1D(:,6),'LineWidth',2)
plot(tsas_1E,xsas_1E(:,6),'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('\Delta{\delta_r} [deg]')
legend('SAS off','SAS 1D','SAS 1E')


figure
subplot(221)
plot(t,y(:,end-1),'LineWidth',2)
hold all
plot(tsas_1D,ysas_1D(:,end-1),'LineWidth',2)
plot(tsas_1E,ysas_1E(:,end-1),'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('\Delta{u_a} [deg]')

subplot(222)
plot(t,y(:,end),'LineWidth',2)
hold all
plot(tsas_1D,ysas_1D(:,end),'LineWidth',2)
plot(tsas_1E,ysas_1E(:,end),'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('\Delta{u_r} [deg]')

subplot(223)
plot(t,x(:,5),'LineWidth',2)
hold all
plot(tsas_1D,xsas_1D(:,5),'LineWidth',2)
plot(tsas_1E,xsas_1E(:,5),'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('\Delta{\delta_a} [deg]')

subplot(224)
plot(t,x(:,6),'LineWidth',2)
hold all
plot(tsas_1D,xsas_1D(:,6),'LineWidth',2)
plot(tsas_1E,xsas_1E(:,6),'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('\Delta{\delta_r} [deg]')
legend('SAS off','SAS 1D','SAS 1E')


figure
subplot(121)
plot(t,x(:,4),'LineWidth',2)
hold all
plot(tsas_1D,xsas_1D(:,4),'LineWidth',2)
plot(tsas_1E,xsas_1E(:,4),'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('r [deg/s]')

subplot(122)
plot(t,y(:,1),'LineWidth',2)
hold all
plot(tsas_1D,ysas_1D(:,1),'LineWidth',2)
plot(tsas_1E,ysas_1E(:,1),'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('r_w [deg/s]')
legend('SAS off','SAS 1D','SAS 1E')
