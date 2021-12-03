%% IDP for Project Mechanism
% Bailey Smoorenburg, Connor McCarthy, Gavin Sheng, Jill Bohnet, Patrick Herke
clear; clc;

w2 = 2*pi;
P4 = -10;

%function ax = frictionIDP(w2, P4)
    tol = 1e-8;
    mu = 0.2;
    R = 0.2; % in
    rw = 1.7/2 % in % half the width of the slider
    rh = 1.0/2 % in % half the height of the slider
    
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
    
    % additional inputs for friction
    f3 = kinData.f3;
    f5 = kinData.f5;

    g = 32.2; %ft/s^2

    % m2 = 0; %slug
    % Ig_2 = 0; %slug*in^2
    % 
    % m3 = 0; %slug
    % Ig_3 = 0; %slug*in^2
    % 
    % m4 = 0; %slug

    m2 = 0.027/32.2; %slug
    Ig_2 = 0.0088/32.2/12; %slug*ft*in

    m3 = 0.1135/32.2; %slug
    Ig_3 = 0.3786/32.2/12; %slug*ft*in

    m4 = .0812/32.2; %slug

    Rcg3 = 3.347; %inches, pythag thm from pic

    % initializing table to store values
    sz = [361 13];
    varNames = ["F12x","F12y","F23","F23x","F23y","F13","F13x","F13y","F34x","F34y","F14","R7F14","T2"];
    varTypes = repmat("double",1,13);
    forcesIDP = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames',varNames);

     for i = 1:361
         stop = 0;
%         if f4*w2 > 0
%             P4 = -10;
%         else
%             P4 = 0;
%         end
        
        %accelerations
        alpha2 = 0;
        alpha3 = h3(i)*alpha2 + h3p(i)*w2^2; %rad/s^2

        a_g2x = 0;
        a_g2y = 0;

        a_g3x = (fg3x(i)*alpha2 + fg3xp(i)*w2^2)/12; %ft/s^2
        a_g3y = (fg3y(i)*alpha2 + fg3yp(i)*w2^2)/12; %ft/s^2

        a_g3x = (f4(i)*alpha2 + f4p(i)*w2^2)/12; %ft/s^2
        a_g3y = 0;

        a_g4y = 0;
        a_g4x = (f4(i)*alpha2 + f4p(i)*w2^2)/12; %ft/s^2
        
        fric = ones(5,1);
        F12 = 0;
        F34 = 0;
        F23n = 0;
        F13n = 0;
        F14 = 0;
        R7F14 = 0;
        
        % input directional indicator
        Din = sign(w2);
        
        % pin joint directional indicators
        D12 = Din;
        D34 = sign(h3(i))*Din;
        
        % pin in slot directional indicators
        if F23n > 0
            D23 = sign(f3(i) + R*(h3(i) - 1))*Din;
        else
            D23 = sign(f3(i) - R*(h3(i) - 1))*Din;
        end
        
        if F13n > 0
            D13 = sign((f3(i)-f5(i)) + R*h3(i))*Din;
        else
            D13 = sign((f3(i)-f5(i)) - R*h3(i))*Din;
        end
        
        % slider in slot directional indicator
        D14 = sign(-f4(i))*Din;
        
        for j = 1:100
            % T12 = mu*R*|F12|*D12 acting on link 1
            T12 = mu*R*F12*D12; %lbf*in
            % T34 = mu*R*|F34|*D34 acting on link 4
            T34 = mu*R*F34*D34; %lbf*in
            % T13 = mu*R*F13*D13 acting on link 1
            T13 = mu*R*F13n*D13; %lbf*in
            % f13 = mu*F13n acting on link 1
            f13 = mu*abs(F13n)*D13; %lbf
            f13x = f13*cos(th3(i) - pi);
            f13y = f13*sin(th3(i) - pi);
            % T23 = mu*R*F13*D13 acting on link 2
            T23 = mu*R*F23n*D23; %lbf*in
            % f23 = mu*F23n acting on link 2
            f23 = mu*abs(F23n)*D23; %lbf
            f23x = f23*cos(th3(i) - pi);
            f23y = f23*sin(th3(i) - pi);
            % f14 acting on link 4
            N14 = [F14; R7F14]\[1 1; rw -rw]; %lbf
            f14 = mu*(abs(N14(1)) + abs(N14(2)))*D14; %lbf
            % the book swaps the signs on T14 depending on >0 but I don;t
            % think that's necessary
            % it also is used if height above and below the block C pin
            % aren't the same
            if N14(1) > 0
                T14(1) = mu*rh*N14(1)*D14; %lbf*in
            else
                T14(1) = mu*rh*N14(2)*D14; %lbf*in
            end
            
            if N14(2) > 0
                T14(2) = mu*rh*N14(2)*D14; %lbf*in
            else
                T14(2) = mu*rh*N14(2)*D14; %lbf*in
            end
            
            % storing the old friction values
            fricOld = fric;
            fric = [(1/sqrt(1/mu^2 + 1))*F12;
                    (1/sqrt(1/mu^2 + 1))*F34;
                    f13;
                    f23;
                    f14;]
           
            relError = norm(fric - fricOld)/norm(fricOld);
            if relError < tol
                fricStore{i} = fric;
                stop = 1;
                break
            end
            
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
            
            jKinem = [m2*a_g2x;
                    m2*a_g2y;
                    m3*a_g3x;
                    m3*a_g3y;
                    m4*a_g4x;
                    m4*a_g4y;
                    Ig_2*alpha2;
                    Ig_3*alpha3 + m3*Rcg3*(cos(th3(i))*a_g3y - sin(th3(i))*a_g3x + g*cos(th3(i)));
                    0];
            
%             jForce = [f23x;
%                     -m2*g + f23y;
%                     -f13x + -f23x;
%                     -m3*g + -f13y + -f23y;
%                     P4 + f14;
%                     -m4*g;
%                     -T12 + T23;
%                     -T34 + -T13 + -T23;
%                     T34 + T14(1) + T14(2)];

            jForce = [0;
                    -m2*g;
                    0;
                    -m3*g;
                    P4;
                    -m4*g;
                    0;
                    0;
                    0];
            
            jF12 = [0;
                    0;
                    0;
                    0;
                    0;
                    0;
                    -T12;
                    0;
                    0];
                
            jF34 = [0;
                    0;
                    0;
                    0;
                    0;
                    0;
                    0;
                    -T34;
                    T34];
            
            jF13 = [0;
                    0;
                    -f13x;
                    -f13y;
                    0;
                    0;
                    0;
                    -T13;
                    0];

            jF23 = [f23x;
                    f23y;
                    -f23x;
                    -f23y;
                    0;
                    0;
                    T23;
                    -T23;
                    0];
            
            jF14 = [0;
                    0;
                    0;
                    0;
                    f14;
                    0;
                    0;
                    0;
                    T14(1) + T14(2)];
                
            J = jKinem - jForce - jF12 - jF34 - jF13 - jF23 - jF14;
            
            x = A\J;
            
            F12 = sqrt(x(1)^2 + x(2)^2); % F12 = sqrt(F12x^2 + F12y^2)
            F34 = sqrt(x(5)^2 + x(6)^2); % F34 = sqrt(F34x^2 + F34y^2)
            F13n = x(4);
            F23n = x(3);
            F14 = x(7);
            R7F14 = x(8);

%             fricOld = fric;
%             fric = (1/sqrt(1/mu^2 + 1)).*[F12; F34; F13; F23; F14];

            % x = [  1,    2,    3,    4,   5,    6,    7,    8,    9]
            % x = [F12x, F12y, F23n, F13n, F34x, F34y, F14, R7F14, T2]
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
        maxes = (find(forcesIDP.(posColNames{i}) == max(forcesIDP.(posColNames{i}))));
        posMinMax.(posColNames{i}) = [find(forcesIDP.(posColNames{i}) == min(forcesIDP.(posColNames{i}))); maxes(1)];
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
    saveName = sprintf('jointForcesFriction_w2_%0.00f_P4_%0.00f.jpg', w2, P4)
    exportgraphics(ax,saveName)
    
    figure(20)
    plot(kinData.theta2, f4*w2, '-')
%end