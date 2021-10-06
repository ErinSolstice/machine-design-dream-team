clear; clc;

% The values from A139 to A216 of excel sheet are not possible due to
% length of theta3

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

sz = [361 5];
varNames = ["theta2","theta3","R3","R4","R5"];
varTypes = ["double","double","double","double","double"];
posData = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames',varNames)

for j = 1:361;
    th2 = (j-1)*pi/180;
    posData.theta2(j) = th2;
    
    tol = 1e-8;
    xNew = [0; 0; 0; 0;];

    for i = 1:100
        xOld = xNew;

        A = [-R3*sin(th3), cos(th3), 1, 0;
          R3*sin(th3), sin(th3), 0, 0;
          R5*sin(th3), 0, 0, -cos(th3);
         -R5*cos(th3), 0, 0, -sin(th3);];

        b = [-(0 + R4 + R3*cos(th3) - R2*cos(th2));
         -(R1 + 0 + R3*sin(th3) - R2*sin(th2));
         -(R2*cos(th2) - R5*cos(th3) - R6*cos(th6));
         -(R2*sin(th2) - R5*sin(th3) - R6*sin(th6));];

        xNew = A\b;

        th3 = th3 + xNew(1);
        R3 = R3 + xNew(2);
        R4 = R4 + xNew(3);
        R5 = R5 + xNew(4);

        relErr = norm(xNew - xOld) / norm(xOld);

        if relErr < tol
            break
        end
    end
    
    posData.theta3(j) = th3;
    posData.R3(j) = R3;
    posData.R4(j) = R4;
    posData.R5(j) = R5;
    
    fprintf("iterations = %d\nth3 = %f\nR3 = %f\nR4 = %f\nR5 = %f\n\n", i, th3*180/pi, R3, R4, R5)
end

filename = 'positionData.xlsx';
writetable(posData,filename,'Sheet',1,'Range','A1')