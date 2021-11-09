%% Program Setup

clear; clc; % clear previous work

% reading in the data
fullTable = readtable('kinematicData.xlsx');


%% 

% defining and intitializing the table
sz = [1001 9];
varNames = ["t","T2","theta2","w2","a2","Ie","dIe_dt","deg","pos"];
varTypes = repmat("double",1,9);
torData = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames',varNames);

% in lbm and convert to slugs
m2 = 0.027/32.2;
m3 = 0.1135/32.2;
m4 = 0.0812/32.2;

% in lbm*in^2 and convert to slug*in^2
Ig2 = 0.0088/32.2;
Ig3 = 0.3786/32.2;

% in rad/s
w2 = 10*(2*pi)*(1/60);
% in lbf
P4 = 10;

% in rad/s
dw2 = 1;
w2 = 0;

% in rad
theta2 = 0;

% in seconds
dt = 0.001;

% interpolation of KCs
interp = @(col, deg, pos) ((fullTable.(col)(pos+1) - fullTable.(col)(pos))*(deg+1-pos) + fullTable.(col)(pos));

for i=1:1001
    deg = theta2*180/pi;
    while deg > 360
        deg = deg - 360;
    end
    while deg < 0
        deg = deg + 360;
    end
    
    pos = floor(deg+1);
    
    torData.deg(i) = deg;
    torData.pos(i) = pos;
    
    h2 = 1;
    h2p = 0;

    h3 = interp("h3", deg, pos);
    h3p = interp("h3p", deg, pos);
    fg3x = interp("fg3x", deg, pos);
    fg3y = interp("fg3y", deg, pos);
    fg3xp = interp("fg3xp", deg, pos);
    fg3yp = interp("fg3yp", deg, pos);

    fg4x = interp("f4", deg, pos);
    fg4xp = interp("f4p", deg, pos);
    
    f4 = interp("f4", deg, pos);
    f4p = interp("f4p", deg, pos);

    Ie = (Ig2*h2^2) + (m3*(fg3x^2 + fg3y^2) + Ig3*h3^2) + (m4*fg4x^2);
    dIe_dt = 2*(Ig2*h2*h2p) + 2*(m3*(fg3x*fg3xp + fg3y*fg3yp) + Ig3*h3*h3p) + 2*(m4*fg4x*fg4xp);
    
    a2 = (T2 - 0.5*dIe_dt*w2^2 - m3*32.2/12*fg3y + P4*f4)/Ie;
    theta2 = theta2 + w2*dt + 0.5*a2*dt^2;
    w2 = w2 + a2*dt;
    
    torData.t(i) = dt*(i);
    torData.T2(i) = T2;
    torData.theta2(i) = theta2;
    torData.w2(i) = w2;
    torData.a2(i) = a2;
    torData.Ie(i) = Ie;
    torData.dIe_dt(i) = dIe_dt;
end

% writing the data to an excel file
filename = 'torqueStartup.xlsx';
writetable(torData,filename,'Sheet',1,'Range','A1')

tspan = 1:100;

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