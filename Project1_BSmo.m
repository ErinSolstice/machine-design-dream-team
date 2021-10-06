%% Project 1 Plotting
% ME 4133 
% 8/26/21
% Bailey Smoorenburg 

%% Declaring variables from imported data
%{
Unit/dimension notes
theta values in radians
h values are dimensionless
R, f, and f-p values have units of length
%}
 
%Create one 360x13 matrix with all values imported from excel
VarMatrix = readmatrix('Project I Plotting Fall 2020.xlsx', 'Range', 'A:M');

%% Assigning proper columns to variables
% Assigning columns to variable names helps keep organized

% Independent variable for all plots, column one of VarMatrix
theta2 = VarMatrix(3:end,1);

% Assigning Position analysis, VarMatrix columns 2 thru 5
theta3 = VarMatrix(3:end,2);
R3 = VarMatrix(3:end,3);
R4 = VarMatrix(3:end,4); 
R5 = VarMatrix(3:end,5);

% Assigning 1st Order Kinematic Coefficients, VarMatrix columns 6 thru 9
h3 = VarMatrix(3:end,6);
f3 = VarMatrix(3:end,7);
f4 = VarMatrix(3:end,8);
f5 = VarMatrix(3:end,9);

% Assigning 2nd Order Kinematic Coefficients, VarMatrix columns 10 thru 13
h3p = VarMatrix(3:end,10);
f3p = VarMatrix(3:end,11);
f4p = VarMatrix(3:end,12);
f5p = VarMatrix(3:end,13);

%% Graph Validation 
% Using the max/min values and zeros of graphs to verify that 1st order is
% a derivative of Position and 2nd Order is a derivitive of 1st Order

% Finding max, min, and max/min index in matrix of Position and 1st Order
for i = 2:9
    [temp_max, theta2max] = max(VarMatrix(3:end,i)); %returns max value and index
    max_matrix(1,i-1) =  temp_max; %mapping max value to a matrix for accessing later
    max_matrix(2,i-1) = theta2max+2; %mapping max index value to a matrix for accessing later
    [temp_min, theta2min] = min(VarMatrix(3:end,i));
    min_matrix(1,i-1) =  temp_min;
    min_matrix(2,i-1) = theta2min+2;
end

% Using index values of min and max to find corresponding zeros in 1st and 2nd Order 
for i = 1:8
   zeros_max(1,i) = VarMatrix(max_matrix(2,i),i+5); 
   zeros_min(1,i) = VarMatrix(min_matrix(2,i),i+5);
   
   % Displaying a string if the value of the found 'zero' is less than
   % arbitrarily chosen 0.01. Should be displayed 8 times if each zero is
   % truly a zero
   if zeros_max(1,i) < 0.01
       disp('Max verified')
   end
   if zeros_min(1,i) < 0.01
       disp('Min verified')
   end
end


%% Plotting
% Position Analysis
figure(1);
plot(theta2, theta3, 'r', theta2, R3, 'g', theta2, R4, 'b', theta2, R5, 'm', 'LineWidth', 2);
title('Position Analysis');
grid on;
% creating axes labels
xlabel('Theta2 (rad)');
ylabel('Outputs');
%creating legend for each variable plotted
legend('Theta3','R3','R4','R5');
%setting the max value for the x axis manually
xlim([0, 2*pi]);
% setting the x-axis grid to count in radians and creating proper labels in
% radians
set(gca, 'XTick', 0: pi/4: 2*pi);
set(gca,'XTickLabel',{'0','pi/4','pi/2','3*pi/4','pi', '5*pi/4', '3*pi/2', '7*pi/4', '2*pi'});

% figure 2 and 3 follow same logic as figure 1

% 1st Order Kinematic Coefficients
figure(2);
plot(theta2, h3, 'r', theta2, f3, 'g', theta2, f4, 'b', theta2, f5, 'm', 'LineWidth', 2);
title('1st Order Kinematic Coefficients');
grid on;
xlabel('Theta2 (rad)');
ylabel('Outputs');
legend('h3','f3(length)','f4(length)','f5(length)');
xlim([0, 2*pi]);
set(gca, 'XTick', 0: pi/4: 2*pi);
set(gca,'XTickLabel',{'0','pi/4','pi/2','3*pi/4','pi', '5*pi/4', '3*pi/2', '7*pi/4', '2*pi'});

% 2nd Order Kinematic Coefficient
figure(3);
plot(theta2, h3p, 'r', theta2, f3p, 'g', theta2, f4p, 'b', theta2, f5p, 'm', 'LineWidth', 2);
title('2nd Order Kinematic Coefficients');
grid on;
xlabel('Theta2 (rad)');
ylabel('Outputs');
legend('h3p','f3p(length)','f4p(length)','f5p(length)');
xlim([0, 2*pi]);
set(gca, 'XTick', 0: pi/4: 2*pi);
set(gca,'XTickLabel',{'0','pi/4','pi/2','3*pi/4','pi', '5*pi/4', '3*pi/2', '7*pi/4', '2*pi'});
