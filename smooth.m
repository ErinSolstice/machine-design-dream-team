function ran = smooth(dw2,P4)
    %% Program Setup

    % requires input of desired w2 value in radians
    % and the applied load as P4

    % reading in the kinematic data
    fullTable = readtable('kinematicData.xlsx');

    % in lbf
    % P4 = 1.5;
    
    % dw2 = 10*2*pi/60;
    f = 2*pi;
    m = dw2*(f/(2*pi));
    tSpan = 1;

    % defining the function for a smooth acceleration
    accel = @(x) ( m*(1 - cos(f*x)) );
    % also defines what the velocity and displacement should look like
    % given continuos integration of the acceleration function
    vel = @(x) ( m*(x - (1/f)*sin(f*x)) );
    displacement = @(x) ( m*(0.5*x^2 + (1/f^2)*cos(f*x) - 1/f^2) );

    %% Finding Worst Startup Position
    sz = [361 5];
    varNames = {'theta2Init','theta2MaxTor','timeMaxTor','maxTor','sheetName'};
    varTypes = repmat("double",1,5);
    varTypes(5) = "string";
    maxTorStartup = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames',varNames);
    startTor = 0;

    for j = 1:361
        theta2 = j-1; % degrees
        maxTorStartup.theta2Init(j) = theta2*pi/180;
        [theta2MaxTor, timeMaxTor, maxTor, sheetName, torData] = startupSim(fullTable,accel,vel,displacement,theta2,dw2,P4);

        maxTorStartup.theta2MaxTor(j) = theta2MaxTor;
        maxTorStartup.timeMaxTor(j) = timeMaxTor;
        maxTorStartup.maxTor(j) = maxTor;
        maxTorStartup.sheetName(j) = sheetName;

        disp(abs(maxTor) >= abs(startTor))
        if abs(maxTor) >= abs(startTor)
            torDataMax = torData;
            startTor = maxTor;
        end
    end

    maxTorSheetName = sprintf('iTh2=%.4f,w2=%.2f,P4=%.2f',theta2,dw2,P4)
    filename = 'maxTor.xlsx';
    writetable(maxTorStartup,filename,'Sheet',maxTorSheetName,'Range','A1')


    %% Writing the data for worst position to an excel file

    index = find(maxTorStartup.maxTor == startTor);
    disTheta2 = maxTorStartup.theta2Init(index);
    sheetName = maxTorStartup.sheetName(index);

    filename = 'torqueStartup.xlsx';
    writetable(torDataMax,filename,'Sheet',sheetName,'Range','A1')

    %% Graphing data for startup torque

    %sheetName = 'initTh2 = 48.0000';
    graphStartup(sheetName,maxTorSheetName,index,dw2,P4)
    
    ran = true;


    %% Startup Simulation Function
    function [theta2MaxTor, timeMaxTor, maxTor, sheetName, torData] = startupSim(fullTable,accel,vel,displacement,theta2,dw2,P4) 
        sheetName = sprintf('iTh2=%.4f,w2=%.2f,P4=%.2f',theta2,dw2,P4)
        theta2 = theta2*pi/180; % convert to radians

        % defining and intitializing the table
        sz = [1501 11];
        varNames = {'t','T2','theta2Sim','theta2Fun','w2Sim','w2Fun','a2','Ie','dIe_dt','deg','pos'};
        varTypes = repmat("double",1,11);
        torData = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames',varNames);

        % in lbm and convert to slugs
        m2 = 0.027/32.2;
        m3 = 0.1135/32.2;
        m4 = 0.0812/32.2;

        % in lbm*in^2 and convert to slug*in^2
        Ig2 = 0.0088/32.2;
        Ig3 = 0.3786/32.2;

        % in rad/s
        w2 = 0;
        a2 = 0;
        t = 0;

        % in rad
        dth2 = 1*pi/180;

        % interpolation of KCs
        interp = @(col, deg, pos) ((fullTable.(col)(pos+1) - fullTable.(col)(pos))*(deg+1-pos) + fullTable.(col)(pos));

        i = 0;
        while t < 1
            i = i + 1;
            if abs(w2) > 100
                fprintf('Stopped at %d iterations and w2 = %f', i, w2)
                break
            end

            deg = theta2*180/pi;
            while deg >= 360
                deg = mod(deg,360);
            end
            while deg < 0
                deg = 360 - mod(-deg,360);
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

            Ie = (Ig2*h2^2) + (m3*(fg3x^2 + fg3y^2) + Ig3*h3^2) + (m4*fg4x^2);
            dIe_dt = 2*(Ig2*h2*h2p) + 2*(m3*(fg3x*fg3xp + fg3y*fg3yp) + Ig3*h3*h3p) + 2*(m4*fg4x*fg4xp);

            % in seconds
            % dt = fzero(@(t)(displacement(t)-(theta2+1*pi/180)),tSpan) - t;
            a2 = accel(t);
            if a2 > 0
                dt = max(roots([0.5*a2, w2, -dth2]));
                if dt > 0.01
                    dt = 0.01;
                end
            else
                dt = 0.001;
            end

            T2 = Ie*a2 + 0.5*dIe_dt*w2^2 + m3*32.2/12*fg3y - P4*fg4x; 

            theta2 = theta2 + w2*dt + 0.5*a2*dt^2;
            w2 = w2 + a2*dt;
            t = t + dt;

            torData.t(i) = t;
            torData.T2(i) = T2;
            torData.theta2Sim(i) = theta2;
            torData.theta2Fun(i) = displacement(t);
            torData.w2Sim(i) = w2;
            torData.w2Fun(i) = vel(t);
            torData.a2(i) = a2;
            torData.Ie(i) = Ie;
            torData.dIe_dt(i) = dIe_dt;
        end
        torData = torData(1:i,{'t','T2','theta2Sim','theta2Fun','w2Sim','w2Fun','a2','Ie','dIe_dt','deg','pos'});

        sMaxTor = max(torData.T2);
        sMinTor = min(torData.T2);
        if abs(sMaxTor) > abs(sMinTor)
            maxTor = sMaxTor;
        else
            maxTor = sMinTor;
        end

        index = find(torData.T2 == maxTor);
        timeMaxTor = torData.t(index);
        theta2MaxTor = torData.theta2Sim(index);
    end

    %% Graphing Function
    function ran = graphStartup(sheetName,maxTorSheetName,index,dw2,P4)
        torData = readtable('torqueStartup.xlsx','Sheet',sheetName);
        maxTorData = readtable('maxTor.xlsx','Sheet',maxTorSheetName);
        tspan = 1:height(torData);

        figName = sprintf('w2 vs time: w2 = %.4f and P4 = %.2f',dw2,P4);
        figure('Name',figName)
        plot(torData.t(tspan), torData.w2Sim(tspan))
        xlabel('time (s)')
        ylabel('w2 (rad/s)')

        figName = sprintf('theta2 vs time: w2 = %.4f and P4 = %.2f',dw2,P4);
        figure('Name',figName)
        plot(torData.t(tspan), torData.theta2Sim(tspan))
        xlabel('t (s)')
        ylabel('theta2 (rad)')

        figName = sprintf('alpha2 vs time: w2 = %.4f and P4 = %.2f',dw2,P4);
        figure('Name',figName)
        plot(torData.t(tspan), torData.a2(tspan))
        xlabel('t (s)')
        ylabel('a2 (rad/s^2)')

        figName = sprintf('Max Torque vs Initial Theta2: w2 = %.4f and P4 = %.2f',dw2,P4);
        figure('Name',figName)
        plot(maxTorStartup.theta2Init, maxTorStartup.maxTor,'-x','MarkerIndices',index)
        xlabel('theta2 (rad)')
        ylabel('torque (lbf*in)')
        
        ran = 1;
    end
end