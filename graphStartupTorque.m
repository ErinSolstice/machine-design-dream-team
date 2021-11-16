tspan = 1:10;

figure(2)
plot(torData.t(tspan), torData.w2(tspan))
xlabel('time (s)')
ylabel('w2 (rad/s)')

figure(3)
plot(torData.theta2(tspan), torData.w2(tspan))
xlabel('theta2 (rad)')
ylabel('w2 (rad/s)')

figure(4)
plot(torData.t(tspan), torData.theta2(tspan))
xlabel('t (s)')
ylabel('theta2 (rad)')

figure(5)
plot(torData.t(tspan), torData.a2(tspan))
xlabel('t (s)')
ylabel('a2 (rad/s^2)')

figure(6)
plot(torData.theta2, torData.Ie)
figure(7)
plot(torData.theta2, torData.dIe_dt)
figure(8)
plot(torData.deg, torData.Ie)
figure(9)
plot(torData.deg, torData.dIe_dt)

figure(2)