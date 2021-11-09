%% IDP for Project Mechanism
% Bailey Smoorenburg, Connor McCarthy, Gavin Sheng, Jill Bohnet, Patrick Herke

clc;
clear;

KinematicDataMatrix = readmatrix('kinematicData.xlsx', 'Range', 'A:Q');%readtable to ref wirh column names

% Givens - R values in inches
R1 = 3.52;
R2 = 0.82;
R3 = KinematicDataMatrix(:,3);
R5 = KinematicDataMatrix(:,5); 

% Theta values, converted into radians
theta2 = KinematicDataMatrix(:,1);
theta3 = KinematicDataMatrix(:,2);
theta5 = theta3;

%inputs for accelerations
f_g3x = KinematicDataMatrix(:,14);
f_g3y = KinematicDataMatrix(:,15);
fp_g3x = KinematicDataMatrix(:,16);
fp_g3y = KinematicDataMatrix(:,17);
f_g4x = KinematicDataMatrix(:,8);
fp_g4x = KinematicDataMatrix(:,12);
h3 = KinematicDataMatrix(:,6);
h3p = KinematicDataMatrix(:,10);

w2 = 10*(2*pi)*(1/60); % in rads per sec 

%Values from Solidworks
% cx = 1.7; %from solidworks model in inches
% cy = 1;
% dcx = cx/2;
% dcy = cy/2;

g = 32.2*12; %in/s^2

m2 = 0.027/32.2; %slug
Ig_2 = 0.0088/32.2; %slug*in^2

m3 = 0.1135/32.2; %slug
Ig_3 = 0.3786/32.2; %slug*in^2

m4 = .0812/32.2; %slug

Rcg3 = 3.347; %inches, pythag thm from pic
P4 = 10; %lbf

iter = length(KinematicDataMatrix(:,1)');

for i = 1:iter
    %accelerations
    alpha2 = 0;
    alpha3 = h3(i)*alpha2 + h3p(i)*w2^2;
    
    a_g2x = 0;
    a_g2y = 0;

    a_g3x = f_g3x(i)*alpha2 + fp_g3x(i)*w2^2;
    a_g3y = f_g3y(i)*alpha2 + fp_g3y(i)*w2^2;

    a_g3x = f_g4x(i)*alpha2 + fp_g4x(i)*w2^2;
    a_g3y = 0;

    a_g4y = 0;
    a_g4x = f_g4x(i)*alpha2 + fp_g4x(i)*w2^2; % fp_g4x = fp4
    
    % A = 9x9
    A = [1, 0, -cos(theta3(i)-(3*pi/2)),                0,           0,    0, 0, 0, 0;
         0, 1, -sin(theta3(i)-(3*pi/2)),                0,           0,    0, 0, 0, 0;
         0, 0, cos(theta3(i)-(3*pi/2)),     cos(theta3(i)-(3*pi/2)), -1,   0, 0, 0, 0;
         0, 0, sin(theta3(i)-(3*pi/2)),     sin(theta3(i)-(3*pi/2)), 0,   -1, 0, 0, 0;
         0, 0,          0,                              0,           1,    0, 0, 0, 0; %why did you have a 1 in this line for F_13?
         0, 0,          0,                              0,           0,    1, 1, 0, 0;
         0, 0, -R2*sin(theta2(i)-theta3(i)),            0,           0,    0, 0, 0, 1; %I don't think you need that +pi/2
         0, 0,         -R3(i),                      -(R3(i)-R5(i)),  0,    0, 0, 0, 0; %ccw vs cw
         0, 0,          0,                              0,           0,    0, 0, 1, 0];
         
    
    J = [m2*a_g2x;
         m2*a_g2y; %center of gravity is at pin + m2*g;
         m3*a_g3x;
         m3*a_g3y + m3*g;
         m4*a_g4x - P4;
         m4*a_g4y - m4*g;
         Ig_2*alpha2;
         Ig_3*alpha3 + (m3*Rcg3*(cos(theta3(i))*a_g3y - sin(theta3(i))*a_g3x));
         0]; %R1*m4*a_g4x you can't do this without also accounting for the other forces on 4
     
    solution(:,i) = A\J;
end
 
%% Plotting 
F_12x = solution(1,:);
F_12y = solution(2,:);
F_23 = solution(3,:);
F_13 = solution(4,:);
F_34x = solution(5,:);
F_34y = solution(6,:);
F_14 = solution(7,:);
R7F_14 = solution(8,:);
T2 = solution(9,:);

figure(1)
subplot(8,1,1)
plot(theta2, F_12x)
title('F_12x');

subplot(8,1,2)
plot(theta2, F_12y)
title('F_12y');

subplot(8,1,3)
plot(theta2, F_23)
title('F_23');

subplot(8,1,4)
plot(theta2, F_13)
title('F_13');

subplot(8,1,5)
plot(theta2, F_34x)
title('F_34x');

subplot(8,1,6)
plot(theta2, F_34y)
title('F_34y');

subplot(8,1,7)
plot(theta2, F_14)
title('F_14');

subplot(8,1,8)
plot(theta2, R7F_14)
title('R7F_14');

figure(2)
plot(theta2, T2)
title('T2');

% b = 8x1 
% b = [F_12x;
%      F_12y;
%      F_23;
%      F_13;
%      F_34x;
%      F_34y;
%      F_14
%      R7F_14;
%      T2];
 
%{
%     A = [1, 0, -1,                0,                  0,               0,                   0,                      0,                  0;
%          0, 1, 0,                 -1,                 0,               0,                   0,                      0,                  0;
%          0, 0, R2*sin(theta2(i)), -R2*cos(theta2(i)), 0,               0,                   0,                      0,                  0;
%          0, 0, 1,                 0,                  1,               0,                   -1,                     0,                  0;
%          0, 0, 0,                 1,                  0,               1,                   0,                      -1,                 0;
%          0, 0, 0,                 0,          -R5(i)*sin(theta5(i)), R5(i)*cos(theta5(i)), R3(i)*sin(theta3(i)), -R3(i)*cos(theta3(i)), 0;
%          0, 0, 0,                 0,                  0,               0,                   1,                      0,                  0;
%          0, 0, 0,                 0,                  0,               0,                   0,                      1,                  1;
%          0, 0, 0,                 0,                  0,               0,                   -dcy,                   dcx,                dcx];

last line of j 
%        m4*(dcx*a_g4y - dcy*a_g4x) + dcy*P4];

%}