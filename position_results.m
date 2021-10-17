clear; clc;

% The values from A139 to A216 of excel sheet are not possible due to length of theta3

% defining known lengths and angles
R1 = 3.52;
th1 = 90*pi/180;
R2 = 0.82;
th4 = 0*pi/180;
R6 = 1.85;
th6 = 63*pi/180;
RCG3 = 1.30;

% Remaining unkown lengths and angles are:
th2 = 0*pi/180;
th3 = 269.3*pi/180;
R3 = 4.07;
R4 = 0.87;
R5 = 1.65;
% th5 = th3;

% defining and intitializing the table
sz = [361 13];
varNames = ["theta2","theta3","R3","R4","R5", "h3","f3","f4","f5", "h3'","f3'","f4'","f5'"];
varTypes = repmat("double",1,13);
posData = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames',varNames)

% generating the data for the position plots
% the for loop uses the previous position data as a guess for the next part
for j = 1:361;
    % calculating theta2 in radians and adding it to the table
    th2 = (j-1)*pi/180;
    posData.theta2(j) = th2;
    
    % setting the tolerance
    tol = 1e-8;
    % setting the correction factor to zero to start
    xNew = [0; 0; 0; 0;];

    % implementing newton rasphon up to 100 times
    for i = 1:100
        xOld = xNew;

        % defining the jacobian
        J = [-R3*sin(th3), cos(th3), 1, 0;
          R3*sin(th3), sin(th3), 0, 0;
          R5*sin(th3), 0, 0, -cos(th3);
         -R5*cos(th3), 0, 0, -sin(th3);];

        % defining the residual
        b = [-(0 + R4 + R3*cos(th3) - R2*cos(th2));
         -(R1 + 0 + R3*sin(th3) - R2*sin(th2));
         -(R2*cos(th2) - R5*cos(th3) - R6*cos(th6));
         -(R2*sin(th2) - R5*sin(th3) - R6*sin(th6));];

        % calculating the correction factor
        xNew = J\b;

        % updating the position variables
        th3 = th3 + xNew(1);
        R3 = R3 + xNew(2);
        R4 = R4 + xNew(3);
        R5 = R5 + xNew(4);

        % calculating the error
        relErr = norm(xNew - xOld) / norm(xOld);

        % ending the loop if the error is below the tolerance
        if relErr < tol
            break
        end
    end
    
    % adding the position data to the table
    posData.theta3(j) = th3;
    posData.R3(j) = R3;
    posData.R4(j) = R4;
    posData.R5(j) = R5;
end

% writing the data to an excel file
filename = 'kinematicData.xlsx';
writetable(posData,filename,'Sheet',1,'Range','A1')