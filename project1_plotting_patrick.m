%% Program Setup

% author: Patrick Herke
% date: 08/24/2021
% ME 4133

clear; clc; % clear previous work

% reading in the data
fullTable = readtable('kinematic_results.xlsx');

% defining colors for the graphs
graphColors = {'#BD3131','#CAC006','#3C7FE6','#40A72A'};
rowMinMax = {'min','max'};


%%  Position Graph
% 

% defining columns for the position graph data
posColNames = {'theta3','R3','R4','R5'};
% initializing the table for the minimums and maximums
posMinMax = table([0;0],[0;0],[0;0],[0;0],'VariableNames',posColNames,'RowNames',rowMinMax);
% finding indices of the local minimums and maximums for the position graph
for i=1:4
    posMinMax.(posColNames{i}) = [find(fullTable.(posColNames{i}) == min(fullTable.(posColNames{i}))); find(fullTable.(posColNames{i}) == max(fullTable.(posColNames{i})))];
end

% defining the figure
figure('Name','Position','position',[10,10,1200,1000])
% plotting theta3, R3, R4, R5 versus theta2
for i=1:4
    plot(fullTable.theta2, fullTable.(posColNames{i}),'-x','MarkerIndices',[posMinMax.(posColNames{i}).'],'color',graphColors{i})
    hold on
end
hold off

% adding plot title
title('Position Analysis')
% creating legend for plot
legend('\theta3_{(rad)}','R3_{(in)}','R4_{(in)}','R5_{(in)}')
% labeling the x & y axes
xlabel('\theta2_{(rad)}')
ylabel('Outputs')
% setting xtick values and labels
xticks(0:pi/4:2*pi)
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})

% adding labels to the marks denoting the maximums and minimums
for i=1:2
    for j=1:4
        % getting the x coordinate of the mark
        xPoint = fullTable.theta2(posMinMax{rowMinMax{i},posColNames{j}});
        % getting the y coordinate of the mark
        yPoint = fullTable.(posColNames{j})(posMinMax{rowMinMax{i},posColNames{j}});
        % combining all the information into a sinlge formatted string for the label
        labelString = sprintf('local %s\n%s=%0.4f @ %s2=%0.4f',rowMinMax{i},posColNames{j},yPoint,char(952),xPoint);
        % adding the text to the graph
        text(xPoint, yPoint, labelString,'VerticalAlignment','top','HorizontalAlignment','center')
    end
end




%% 1st Order Kinematic Coefficients Graph
% 

% defining columns for the First Order Kinematic Coefficients graph data
firstOrderColNames = {'h3','f3','f4','f5'};
% initializing the table for the minimums and maximums
firstOrderMinMax = table([0;0],[0;0],[0;0],[0;0],'VariableNames',firstOrderColNames,'RowNames',{'min','max'});
% finding indices of the local minimums and maximums for the First Order Kinematic Coefficients graph
for i=1:4
    firstOrderMinMax.(firstOrderColNames{i}) = [find(fullTable.(firstOrderColNames{i}) == min(fullTable.(firstOrderColNames{i}))); find(fullTable.(firstOrderColNames{i}) == max(fullTable.(firstOrderColNames{i})))];
end

% defining figure
figure('Name','1st Order','position',[10,10,1200,1000])
% plotting theta3, R3, R4, R5 versus theta2
% plotting the mins and maxes
for i=1:4
    plot(fullTable.theta2, fullTable.(firstOrderColNames{i}),'-x','MarkerIndices',[firstOrderMinMax.(firstOrderColNames{i}).'],'color',graphColors{i})
    hold on
end
% plotting the zeros
for i=1:4
    plot(fullTable.theta2(posMinMax.(posColNames{i})),[0;0],'x','color',graphColors{i})
    hold on
end
hold off

% adding plot title
title('1st Order Kinematic Coefficients')
% creating legend for plot
legend('h3_{(-)}','f3_{(length)}','f4_{(length)}','f5_{(length)}')
% labeling the x & y axes
xlabel('\theta2_{(rad)}')
ylabel('Outputs')
% setting xtick values and labels
xticks(0:pi/4:2*pi)
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})

% adding labels to the marks denoting the maximums and minimums
% setting the alignments of the labels relative to the marks
firstVertAlign = {'top','top','top','bottom';'top','top','top','top'};
for i=1:2
    for j=1:4
        % getting the x coordinate of the mark
        xPoint = fullTable.theta2(firstOrderMinMax{rowMinMax{i},firstOrderColNames{j}});
        % getting the y coordinate of the mark
        yPoint = fullTable.(firstOrderColNames{j})(firstOrderMinMax{rowMinMax{i},firstOrderColNames{j}});
        % combining all the information into a sinlge formatted string for the label
        labelString = sprintf('local %s\n%s=%0.4f @ %s2=%0.4f\n',rowMinMax{i},firstOrderColNames{j},yPoint,char(952),xPoint);
        % adding the text to the graph
        text(xPoint, yPoint, labelString,'VerticalAlignment',firstVertAlign{i,j},'HorizontalAlignment','center')
    end
end

% adding labels to the marks denoting the zeros
% setting the alignments of the labels relative to the marks
firstVertAlign = {'bottom','top','top','bottom';'bottom','bottom','top','top'};
firstHorizAlign = {'center','center','center','center';'center','left','center','center'};
for i=1:2
    for j=1:4
        % getting the x coordinate of the mark
        xPoint = fullTable.theta2(posMinMax{rowMinMax{i},posColNames{j}});
        % setting the y to 0
        yPoint = 0;
        % combining all the information into a sinlge formatted string for the label
        labelString = sprintf('\n%s=%0.1f @ %s2=%0.4f\n',firstOrderColNames{j},yPoint,char(952),xPoint);
        % adding the text to the graph
        text(xPoint, yPoint, labelString,'VerticalAlignment',firstVertAlign{i,j},'HorizontalAlignment',firstHorizAlign{i,j})
    end
end




%% 2nd Order Kinematic Coefficients Graph
% 

% defining columns for the Second Order Kinematic Coefficients graph data
secondOrderColNames = {'h3p','f3p','f4p','f5p'};
% initializing the table for the minimums and maximums
secondOrderMinMax = table([0;0],[0;0],[0;0],[0;0],'VariableNames',secondOrderColNames,'RowNames',{'min','max'});
% finding indices of the local minimums and maximums for the secondOrderition graph
for i=1:4
    secondOrderMinMax.(secondOrderColNames{i}) = [find(fullTable.(secondOrderColNames{i}) == min(fullTable.(secondOrderColNames{i}))); find(fullTable.(secondOrderColNames{i}) == max(fullTable.(secondOrderColNames{i})))];
end

% defining figure
figure('Name','2nd Order','position',[10,10,1200,1000])
% plotting theta3, R3, R4, R5 versus theta2
for i=1:4
    plot(fullTable.theta2, fullTable.(secondOrderColNames{i}),'color',graphColors{i})
    hold on
end
% plotting the zeros
for i=1:4
    plot(fullTable.theta2(firstOrderMinMax.(firstOrderColNames{i})),[0;0],'x','color',graphColors{i})
    hold on
end
hold off

% adding plot title
title('2nd Order Kinematic Coefficients')
% creating legend for plot
legend("h3'_{(-)}","f3'_{(length)}","f4'_{(length)}","f5'_{(length)}")
% labeling the x & y axes
xlabel('\theta2_{(rad)}')
ylabel('Outputs')
% setting xtick values and labels
xticks(0:pi/4:2*pi)
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})

% adding labels to the marks denoting the maximums and minimums
% setting the alignments of the labels relative to the marks
secondVertAlign = {'bottom','top','bottom','bottom';'bottom','bottom','top','top'};
secondHorizAlign = {'right','center','center','center';'center','left','center','left'};
for i=1:2
    for j=1:4
        % getting the x coordinate of the mark
        xPoint = fullTable.theta2(firstOrderMinMax{rowMinMax{i},firstOrderColNames{j}});
        % setting the y to 0
        yPoint = 0;
        % combining all the information into a sinlge formatted string for the label
        labelString = sprintf('\n%s=%0.1f @ %s2=%0.4f\n',secondOrderColNames{j},yPoint,char(952),xPoint);
        % adding the text to the graph
        text(xPoint, yPoint, labelString,'VerticalAlignment',secondVertAlign{i,j},'HorizontalAlignment',secondHorizAlign{i,j})
    end
end