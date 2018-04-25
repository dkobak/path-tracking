clear all
load preprocessedData.mat

for sub = 1:size(trainingData,1)
    for day = 1:5
        for tr = 1:length(trainingData{sub,day})
            cursor = trainingData{sub,day}{tr}.cursorPosition;
            borders = trainingData{sub,day}{tr}.pathBorders;
            accuracy(sub,day,tr) = sum(cursor<borders(:,1) & cursor>borders(:,2)) / length(cursor) * 100;
        end
    end
end

expertColor = [55,126,184]/256;

figure('Position', [100 100 800 300])
subplot(121)
hold on
plot(mean(accuracy,3)', 'Color', expertColor)
plot(mean(mean(accuracy,3),1),'.-','LineWidth',3,'color', [0 0 0],'MarkerSize',20)
xlabel('Training day')
ylabel('Time inside the path (%)')
title('Overall performance')

rescaled = mean(accuracy,3) - mean(accuracy(:,1,:),3);

subplot(122)
hold on
plot(rescaled', 'Color', expertColor)
plot(mean(rescaled,1),'.-','LineWidth',3,'color', [0 0 0],'MarkerSize',20)
xlabel('Training day')
title('Performance relative to day 1')

day1 = mean(accuracy(:,1,:),3);
day5 = mean(accuracy(:,5,:),3);
[p,~,stats] = signrank(day1, day5);
d = (mean(day5)-mean(day1))/std([day1-mean(day1); day5-mean(day5)]);
disp(' ')
display(['Day 1: ' num2str(mean(day1)) ' +- ' num2str(std(day1))])
display(['Day 5: ' num2str(mean(day5)) ' +- ' num2str(std(day5))])
display(['Signed-rank p=' num2str(p) ', z=', num2str(abs(stats.zval))])
display(['Cohen`s d=' num2str(d)])
