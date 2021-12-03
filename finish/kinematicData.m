function posData = kinematicData()
    clear; clc;

    % The values from A139 to A216 of excel sheet are not possible due to length of theta3

    % defining known lengths and angles
    R1 = 3.52;
    th1 = 90*pi/180;
    R2 = 0.82;
    th4 = 0*pi/180;
    R6 = 1.85;
    th6 = 63*pi/180;
    RCG3 = 2.28;

    % Remaining unkown lengths and angles are:
    th2 = 0*pi/180;
    th3 = 269.3*pi/180;
    R3 = 4.07;
    R4 = 0.87;
    R5 = 1.65;
    % th5 = th3;

    % defining and intitializing the table
    sz = [361 17];
    varNames = ["theta2","theta3","R3","R4","R5", "h3","f3","f4","f5", "h3p","f3p","f4p","f5p", "fg3x","fg3y","fg3xp","fg3yp"];
    varTypes = repmat("double",1,17);
    posData = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames',varNames);

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

    %%First Order Kinematic Coefficients
    for i=1:361
        %%creating variables for the data
        theta2 = posData.theta2(i);
        theta3 = posData.theta3(i);
        R3 = posData.R3(i);
        R4 = posData.R4(i);
        R5 = posData.R5(i);


        %%writing in the Jacobian
        J=[-R3*sin(theta3),cos(theta3),1,0;
           R3*cos(theta3),sin(theta3),0,0;
           R5*sin(theta3),0,0,-cos(theta3);
           -R5*cos(theta3),0,0,-sin(theta3)];

        %%writing in "b"
        b=[-R2*sin(theta2);
            R2*cos(theta2);
            R2*sin(theta2);
            -R2*cos(theta2)];

        %%calculating the first order KC's
        x=J\b;

        posData.h3(i) = x(1,1);
        posData.f3(i) = x(2,1);
        posData.f4(i) = x(3,1);
        posData.f5(i) = x(4,1);
    end

    %%First Order Kinematic Coefficients
    for i=1:361
        %%creating variables for the data
        %%position
        theta2 = posData.theta2(i);
        theta3 = posData.theta3(i);
        R3 = posData.R3(i);
        R4 = posData.R4(i);
        R5 = posData.R5(i);

        %%first order kinemtaics
        h3 = posData.h3(i);
        f3 = posData.f3(i);
        f4 = posData.f4(i);
        f5 = posData.f5(i);

        %%writing in the Jacobian, remains the same from previous
        J = [-R3*sin(theta3),cos(theta3),1,0;
           R3*cos(theta3),sin(theta3),0,0;
           R5*sin(theta3),0,0,-cos(theta3);
           -R5*cos(theta3),0,0,-sin(theta3)];

        %%writing in "b"
        b = [-(R2*cos(theta2) - R3*h3^2*cos(theta3) - 2*f3*h3*sin(theta3));
            -(R2*sin(theta2) - R3*h3^2*sin(theta3) + 2*f3*h3*cos(theta3));
            -(-R2*cos(theta2) + R5*h3^2*cos(theta3) + 2*f5*h3*sin(theta3));
            -(-R2*sin(theta2) + R5*h3^2*sin(theta3) - 2*f5*h3*cos(theta3))];

        %%calculating the first order KC's
        x = J\b;

        posData.h3p(i) = x(1,1);
        posData.f3p(i) = x(2,1);
        posData.f4p(i) = x(3,1);
        posData.f5p(i) = x(4,1);
    end

    for i=1:361
        %%creating variables for the data
        %%position
        theta2 = posData.theta2(i);
        theta3 = posData.theta3(i);
        RCG3 = 2.28;

        h3 = posData.h3(i);
        f4 = posData.f4(i);

        h3p = posData.h3p(i);
        f4p = posData.f4p(i);

        fg3x = f4 - RCG3*h3*sin(theta3);
        fg3y = RCG3*h3*cos(theta3);

        fg3xp = f4p + RCG3*(-h3p*sin(theta3) - h3^2*cos(theta3));
        fg3yp = RCG3*(h3p*cos(theta3) - h3^2*sin(theta3));

        posData.fg3x(i) = fg3x;
        posData.fg3y(i) = fg3y;
        posData.fg3xp(i) = fg3xp;
        posData.fg3yp(i) = fg3yp;
    end

    % writing the data to an excel file
    filename = 'kinematicData.xlsx';
    writetable(posData,filename,'Sheet',1,'Range','A1')
end