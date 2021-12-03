function torqueConstant(w2,P4)
    %% Program Setup

    % reading in the data
    fullTable = readtable('kinematicData.xlsx');


    %% 

    % defining and intitializing the table
    sz = [361 2];
    varNames = ["theta2","T2"];
    varTypes = repmat("double",1,2);
    torDataConst = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames',varNames);

    % in lbm and convert to slugs
    m2 = 0.027/32.2;
    m3 = 0.1135/32.2;
    m4 = 0.0812/32.2;

    % in lbm*in^2 and convert to slug*ft*in
    Ig2 = 0.0088/32.2/12;
    Ig3 = 0.3786/32.2/12;

    g = 32.2 % ft/s^2

    for i = 1:361
        h2 = 1;
        h2p = 0;

        h3 = fullTable.h3(i);
        h3p = fullTable.h3p(i);
        fg3x = fullTable.fg3x(i);
        fg3y = fullTable.fg3y(i);
        fg3xp = fullTable.fg3xp(i);
        fg3yp = fullTable.fg3yp(i);

        fg4x = fullTable.f4(i);
        fg4xp = fullTable.f4p(i);

    %     if fg4x*w2 > 0
    %         P4 = -1.5;
    %     else
    %         P4 = 0;
    %     end

        dIe_dt = 2*(Ig2*h2*h2p) + 2*(m3*(fg3x*fg3xp + fg3y*fg3yp)/12 + Ig3*h3*h3p) + 2*(m4*fg4x*fg4xp)/12;

        T2 = 0.5*dIe_dt*w2^2 + m3*g*fg3y - P4*fg4x;

        torDataConst.theta2(i) = (i-1)*pi/180;
        torDataConst.T2(i) = T2;
    end

    figure(3)
    plot(torDataConst.theta2,torDataConst.T2)
    xlabel('theta2 (rad)')
    ylabel('torque (lbs)')
    legend('torque')

    % writing the data to an excel file
    filename = 'torqueConstant.xlsx';
    writetable(torDataConst,filename,'Sheet',1,'Range','A1')
end