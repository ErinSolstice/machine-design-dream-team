%% IDP for Project Mechanism
% Bailey Smoorenburg, Connor McCarthy, Gavin Sheng, Jill Bohnet, Patrick Herke

function ax = frictionIDP(w2, P4)
    tol = 1e-8;
    mu = 0.2;
    R = 0.1; % in
    
    kinData = readtable('kinematicData.xlsx');

    % Givens - R values in inches
    R1 = 3.52;
    R2 = 0.82;
    R3 = kinData.R3;
    R5 = kinData.R5; 

    % Theta values, converted into radians
    th2 = kinData.theta2;
    th3 = kinData.theta3;

    %inputs for accelerations
    fg3x = kinData.fg3x;
    fg3y = kinData.fg3y;
    fg3xp = kinData.fg3xp;
    fg3yp = kinData.fg3yp;
    f4 = kinData.f4;
    f4p = kinData.f4p;
    h3 = kinData.h3;
    h3p = kinData.h3p;

    g = 32.2; %in/s^2

    % m2 = 0; %slug
    % Ig_2 = 0; %slug*in^2
    % 
    % m3 = 0; %slug
    % Ig_3 = 0; %slug*in^2
    % 
    % m4 = 0; %slug

    m2 = 0.027/32.2; %slug
    Ig_2 = 0.0088/32.2; %slug*in^2

    m3 = 0.1135/32.2; %slug
    Ig_3 = 0.3786/32.2; %slug*in^2

    m4 = .0812/32.2; %slug

    Rcg3 = 3.347; %inches, pythag thm from pic

    sz = [361 13];
    varNames = ["F12x","F12y","F23","F23x","F23y","F13","F13x","F13y","F34x","F34y","F14","R7F14","T2"];
    varTypes = repmat("double",1,13);
    forcesIDP = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames',varNames);

    for i = 1:361
        %accelerations
        alpha2 = 0;
        alpha3 = h3(i)*alpha2 + h3p(i)*w2^2;

        a_g2x = 0;
        a_g2y = 0;

        a_g3x = (fg3x(i)*alpha2 + fg3xp(i)*w2^2)/12;
        a_g3y = (fg3y(i)*alpha2 + fg3yp(i)*w2^2)/12;

        a_g3x = (f4(i)*alpha2 + f4p(i)*w2^2)/12;
        a_g3y = 0;

        a_g4y = 0;
        a_g4x = (f4(i)*alpha2 + f4p(i)*w2^2)/12; % fp_g4x = fp4
        
        fric = zeros(1,1);
        FMag = zeros(2,1);
        F12 = 0;
        F34 = 0;
        F23 = 0;
        
        % pin joint directional indicators
        D12 = sign(w2);
        D34 = sign(h3(i))*sign(w2);
        
        % pin in slot directional indicators
        if F23 > 0
            D23 = 0;
        else
            D23 = 0;
        end
        
        D13 = 0;
        
        % slider in slot directional indicator
        D14 = 0;
        
        for j = 1:100
            % A = 9x9
            A = [1, 0, -cos(th3(i)-(3*pi/2)),                0,           0,    0, 0, 0, 0;
                 0, 1, -sin(th3(i)-(3*pi/2)),                0,           0,    0, 0, 0, 0;
                 0, 0, cos(th3(i)-(3*pi/2)),     cos(th3(i)-(3*pi/2)), -1,   0, 0, 0, 0;
                 0, 0, sin(th3(i)-(3*pi/2)),     sin(th3(i)-(3*pi/2)), 0,   -1, 0, 0, 0;
                 0, 0,          0,                              0,           1,    0, 0, 0, 0;
                 0, 0,          0,                              0,           0,    1, 1, 0, 0;
                 0, 0, -R2*cos(th3(i)-th2(i)),     0,           0,    0, 0, 0, 1;
                 0, 0,         R3(i),                      R3(i)-R5(i),  0,    0, 0, 0, 0;
                 0, 0,          0,                              0,           0,    0, 0, 1, 0];


            % T12 = mu*R*|F12|*D12
            T12 = mu*R*F12*D12;
            % T34 = mu*R*|F34|*D34;
            T34 = mu*R*F34*D34;
            
            J = [m2*a_g2x;
                 m2*a_g2y + m2*g; %center of gravity is at pin + m2*g;
                 m3*a_g3x;
                 m3*a_g3y + m3*g;
                 m4*a_g4x - P4;
                 m4*a_g4y + m4*g;
                 Ig_2*alpha2/12 + T12;
                 Ig_3*alpha3/12 + m3*Rcg3*(cos(th3(i))*a_g3y - sin(th3(i))*a_g3x + m3*g*cos(th3(i))) - T34;
                 0 + T34]; %change you can't do this without also accounting for the other forces on 4
            
            x = A\J;
            
            F12 = sqrt(x(1)^2 + x(2)^2); % F12 = sqrt(F12x^2 + F12y^2)
            F34 = sqrt(x(5)^2 + x(6)^2); % F34 = sqrt(F34x^2 + F34y^2)
            
            fricOld = fric;
            fric = (1/sqrt(1/mu^2 + 1)).*[F12; F34;]

            
%             fricOld = fric;
%             fric = (1/sqrt(1/mu^2 + 1)).*[F12; F34; F13; F23; F14];

            % x = [  1,    2,   3,   4,   5,    6,    7,    8,    9]
            % x = [F12x, F12y, F23, F13, F34x, F34y, F14, R7F14, T2]
            
            relError = norm(fric - fricOld)/norm(fricOld);
            if relError < tol
                break
            end
        end
        disp(j)

        forcesIDP.F12x(i) = x(1);
        forcesIDP.F12y(i) = x(2);

        forcesIDP.F23(i) = x(3);
        forcesIDP.F23x(i) = x(3)*cos(th3(i)-(3*pi/2));
        forcesIDP.F23y(i) = x(3)*sin(th3(i)-(3*pi/2));

        forcesIDP.F13(i) = x(4);
        forcesIDP.F13x(i) = x(4)*cos(th3(i)-(3*pi/2));
        forcesIDP.F13y(i) = x(4)*sin(th3(i)-(3*pi/2));

        forcesIDP.F34x(i) = x(5);
        forcesIDP.F34y(i) = x(6);

        forcesIDP.F14(i) = x(7);
        forcesIDP.R7F14(i) = x(8);

        forcesIDP.T2(i) = x(9);

    end
    
    filename = 'forcesIDP.xlsx';
    writetable(forcesIDP,filename,'Sheet',1,'Range','A1')

    %% Plotting

    % defining colors for the graphs
    graphColors = {'#B58900','#cb4b16','#dc322f','#d33682','#6c71c4','#268bd2','#2aa198','#859900'};
    rowMinMax = {'min','max'};


    %%  Position Graph

    % defining columns for the position graph data

    % initializing the table for the minimums and maximums
    sz = [2 8];
    posColNames = {'F12x','F12y','F23','F13','F34x','F34y','F14','T2'};
    varTypes = repmat("double",1,8);
    posMinMax = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames',posColNames, 'RowNames',rowMinMax);

    % finding indices of the local minimums and maximums for the position graph
    for i=1:8
        posMinMax.(posColNames{i}) = [find(forcesIDP.(posColNames{i}) == min(forcesIDP.(posColNames{i}))); find(forcesIDP.(posColNames{i}) == max(forcesIDP.(posColNames{i})))];
    end

    % defining the figure
    figure('Name','Position','position',[10,10,1200,1000])
    % plotting theta3, R3, R4, R5 versus theta2
    for i=1:8
        plot(kinData.theta2, forcesIDP.(posColNames{i}),'-x','MarkerIndices',[posMinMax.(posColNames{i}).'],'color',graphColors{i})
        hold on
    end
    hold off

    % adding plot title
    tit = sprintf('Joint Force Analysis at w2 = %0.00f and P4 = %0.00f', w2, P4)
    title(tit)
    % creating legend for plot
    legend(posColNames)
    % labeling the x & y axes
    xlabel('\theta2_{(rad)}')
    ylabel('Joint Forces')
    % setting xtick values and labels
    xticks(0:pi/4:2*pi)
    xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
    % adding gridlines
    grid on
    grid minor
    % setting the limits of the axis
    xlim([0,2*pi])

    % adding labels to the marks denoting the maximums and minimums
    % for i=1:2
    %     for j=1:8
    %         % getting the x coordinate of the mark
    %         xPoint = kinData.theta2(posMinMax{rowMinMax{i},posColNames{j}});
    %         % getting the y coordinate of the mark
    %         yPoint = kinData.(posColNames{j})(posMinMax{rowMinMax{i},posColNames{j}});
    %         % combining all the information into a sinlge formatted string for the label
    %         labelString = sprintf('local %s\n%s=%0.4f @ %s2=%0.4f',rowMinMax{i},posColNames{j},yPoint,char(952),xPoint);
    %         % adding the text to the graph
    %         text(xPoint, yPoint, labelString,'VerticalAlignment','top','HorizontalAlignment','center')
    %     end
    % end

    % saving the graph
    ax = gca;
    saveName = sprintf('jointForces_w2_%0.00f_P4_%0.00f.jpg', w2, P4)
    exportgraphics(ax,saveName)
end