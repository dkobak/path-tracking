clear all
load 'preprocessedData.mat'

n = zeros(length(expData), 10);
accuracy = zeros(3, length(expData), 10);

% compute accuracy for each subject and trial
for subject = 1:length(expData)
    for trial = 1:length(expData{1})
        searchlight = round(expData{subject}{trial}.searchlightLength/10);
        n(subject, searchlight) = n(subject, searchlight) + 1;
        
        cursor = expData{subject}{trial}.cursorPosition;
        borders = expData{subject}{trial}.pathBorders;
        midline = expData{subject}{trial}.pathMidline;
        
        accuracy(n(subject, searchlight), subject, searchlight) =  ...
            sum(cursor<borders(:,1) & cursor>borders(:,2)) / length(cursor) * 100;
    end
end
save('accuracy.mat', 'accuracy')

% Compare naive and expert on full searchlight
naive  = squeeze(mean(accuracy(:,naiveGroup,10),1));
expert = squeeze(mean(accuracy(:,expertGroup,10),1));
[p,~,stats] = ranksum(naive, expert);
disp(' ')
d = abs(mean(naive)-mean(expert))/std([naive-mean(naive) expert-mean(expert)]);
display(['Expert s=100: ' num2str(mean(expert)) ' +- ' num2str(std(expert))])
display(['Naive  s=100: ' num2str(mean(naive)) ' +- ' num2str(std(naive))])
display(['Ranksum p=' num2str(p) ', z=', num2str(abs(stats.zval))])
display(['Cohen`s d=' num2str(d)])

% Compare naive and expert on the smallest searchlight
naive  = squeeze(mean(accuracy(:,naiveGroup,1),1));
expert = squeeze(mean(accuracy(:,expertGroup,1),1));
[p,~,stats] = ranksum(naive, expert);
d = abs(mean(naive)-mean(expert))/std([naive-mean(naive) expert-mean(expert)]);
disp(' ')
display(['Expert s=10: ' num2str(mean(expert)) ' +- ' num2str(std(expert))])
display(['Naive  s=10: ' num2str(mean(naive)) ' +- ' num2str(std(naive))])
display(['Ranksum p=' num2str(p) ', z=', num2str(abs(stats.zval))])
display(['Cohen`s d=' num2str(d)])

% compute horizons for each subject
for subj = 1:length(expData)
    y = squeeze(accuracy(:,subj,:));
    x = repmat(10:10:100, [3 1]);

    beta(subj,:) = changepoint(y, 10:10:100);
    
    expFun = @(b,x)(b(2)*exp(x.*b(3))+b(1));
    expBeta(subj,:) = nlinfit(x(:), y(:), expFun, [mean(mean(y(:,8:10))) -50 -1/5]);    
end

asymptote = squeeze(mean(mean(accuracy(:,:,8:10),3),1));
slope = beta(:,1) * 100 / 30;
horizonChangepoint = beta(:,3);
horizonChangepoint = horizonChangepoint / 100 * 30;
interceptChangepoint = beta(:,2);
% 5 is an arbitrary number. It's like half-life, but only until 5-fold decrease
% adding 10 because that's the shortest searchlight
horizonExp = -1./expBeta(:,3) * log(5) + 10;
horizonExp = horizonExp / 100 * 30;

% % Plot raw data for all subjects with fits
% [~, subjSort] = sort(asymptote, 'descend');
% figure('Position', [100 100 1600 800])
% x = 10:10:100;
% xx = 10:1:100;
% for subj = 1:size(accuracy, 2)
%     subplot(8,8,subj)
%     if expertGroup(subjSort(subj))
%         title(['expert ' num2str(subjSort(subj))])
%     else
%         title(['naive ' num2str(subjSort(subj))])
%     end
%     hold on
%     
%     y = squeeze(accuracy(:,subjSort(subj),:));
%     plot(x, y', 'b.')
%         
%     ylim([0 100])
%     xlim([0 110])
%     xticks([])
%     yticks([])
%     
%     b = beta(subjSort(subj),:);
%     yhat = b(2)+b(1)*xx;
%     yhat(xx>b(3)) = b(2)+b(1)*b(3);
%     plot(xx, yhat, 'k')
%     plot(b(3)*[1 1], [0 b(2)+b(1)*b(3)], 'k--')
% end

naiveColor = [228,26,28]/256;
expertColor = [55,126,184]/256;
cm = (10:10:100) / 100 * 30;

% Actual figure
figure('Position', [100 100 1200 700])
subplot(232)
hold on
plot(cm, squeeze(mean(accuracy(:,expertGroup,:), 1))', 'Color', expertColor)
plot(cm, squeeze(mean(accuracy(:,naiveGroup,:), 1))', 'Color', naiveColor)
axis([0 30 20 100])
xticks([0 6 12 18 24 30])
yticks(20:10:100)
xlabel('Searchlight length (cm)')
ylabel('Time inside the path (%)')

% Rescaled accuracies
acc_rescaled = squeeze(mean(accuracy, 1));
acc_rescaled = acc_rescaled - acc_rescaled(:,1);
acc_rescaled = acc_rescaled ./ (mean(acc_rescaled(:,8:10),2)-acc_rescaled(:,1));
subplot(233)
hold on
axis([0 30 -.2 1.3])
xticks([0 6 12 18 24 30])
yticks(0:.2:1.2)
xlabel('Searchlight length (cm)')
ylabel('Time inside the path (rescaled)')
plot([0 30], [1 1], 'k--')
plot([0 30], [0 0], 'k--')

naiveGroupReduced = naiveGroup;
naiveGroupReduced([6 17] + sum(expertGroup)) = false;    % excluding two flat subjects
plot(cm, acc_rescaled(expertGroup,:)', 'Color', expertColor * 0.5 + [1 1 1] * 0.5)
plot(cm, acc_rescaled(naiveGroupReduced, :)', 'Color', naiveColor * 0.5 + [1 1 1] * 0.5)

% 95% confidence intervals
mu = mean(acc_rescaled(expertGroup,:),1);
sigma = std(acc_rescaled(expertGroup,:),0,1)/sqrt(sum(expertGroup));
errorbar(cm, mu, 1.96*sigma, '.-', 'Color', expertColor, 'LineWidth', 2, 'MarkerSize', 20)
mu = mean(acc_rescaled(naiveGroupReduced,:),1);
sigma = std(acc_rescaled(naiveGroupReduced,:),0,1)/sqrt(sum(naiveGroupReduced));
errorbar(cm, mu, 1.96*sigma, '.-', 'Color', naiveColor, 'LineWidth', 2, 'MarkerSize', 20)

% Holm-Bonferroni
for i=1:10
    p(i) = ranksum(acc_rescaled(expertGroup,i), acc_rescaled(naiveGroupReduced,i));
end
[pcorr, ind] = sort(p(2:7), 'ascend');
for i=1:length(ind)
    pcorr(i) = pcorr(i) * (length(ind) + 1 - i);
end
p(ind+1) = pcorr;
disp(' ')
display('Ranksum p-values after Holm-Bonferroni:')
display(num2str(p))
text(11.2, 1.2, '**')
text(14.2, 1.2, '**')

% Scatter plots
subplot(234)
hold on
scatter(asymptote(expertGroup), horizonChangepoint(expertGroup), 25, expertColor, 'filled')
scatter(asymptote(naiveGroup), horizonChangepoint(naiveGroup), 25, naiveColor, 'filled')
ylim([0 30])
yticks(0:3:30)
xlim([20 100])
xticks(20:20:100)
xlabel('Time inside the path at the asymptote (%)')
ylabel('Planning horizon, changep. fit (cm)')

% Compare naive and expert horizons
naive  = horizonChangepoint(naiveGroup);
expert = horizonChangepoint(expertGroup);
[p,~,stats] = ranksum(naive, expert);
d = abs(mean(naive)-mean(expert))/std([naive-mean(naive); expert-mean(expert)]);
disp(' ')
display(['Expert horizon: ' num2str(mean(expert)) ' +- ' num2str(std(expert)) '. Median: ' num2str(median(expert))])
display(['Naive  horizon: ' num2str(mean(naive)) ' +- ' num2str(std(naive)) '. Median: ' num2str(median(naive))])
display(['Ranksum p=' num2str(p) ', z=', num2str(abs(stats.zval))])
display(['Cohen`s d=' num2str(d)])
[r, pCorr] = corr(asymptote(:), horizonChangepoint, 'type', 'Spearman');
display(['Spearman corr=' num2str(r) ', p=' num2str(pCorr)])
text(30,20, ['R=' num2str(r,2) ' ***'])
nrep=1000;
rng(42);
for rep=1:nrep
    bootdiff(rep) = median(expert(randi(length(expert), length(expert), 1))) - median(naive(randi(length(naive), length(naive), 1))); 
end
display(['95% CI for the difference in medians: ' num2str(prctile(bootdiff, 2.5)) '--' num2str(prctile(bootdiff, 97.5))])

% [r,pCorr]=corr([asymptote(naiveGroup) - mean(asymptote(naiveGroup)) asymptote(expertGroup) - mean(asymptote(expertGroup))]', [horizonChangepoint(naiveGroup) - mean(horizonChangepoint(naiveGroup)); horizonChangepoint(expertGroup) - mean(horizonChangepoint(expertGroup))], 'type', 'Spearman')
% figure
% scatter([asymptote(naiveGroup) - mean(asymptote(naiveGroup)) asymptote(expertGroup) - mean(asymptote(expertGroup))]', [horizonChangepoint(naiveGroup) - mean(horizonChangepoint(naiveGroup)); horizonChangepoint(expertGroup) - mean(horizonChangepoint(expertGroup))])


% Using exponential fit
subplot(235)
hold on
scatter(asymptote(expertGroup), horizonExp(expertGroup), 25, expertColor, 'filled')
scatter(asymptote(naiveGroup), horizonExp(naiveGroup), 25, naiveColor, 'filled')
ylim([0 30])
yticks(0:3:30)
xlim([20 100])
xticks(20:20:100)
xlabel('Time inside the path at the asymptote (%)')
ylabel('Planning horizon, exp. fit (cm)')

% Compare naive and expert horizons
naive  = horizonExp(naiveGroup);
expert = horizonExp(expertGroup);
[p,~,stats] = ranksum(naive, expert);
d = abs(mean(naive)-mean(expert))/std([naive-mean(naive); expert-mean(expert)]);
disp(' ')
display(['Expert horizon exp: ' num2str(mean(expert)) ' +- ' num2str(std(expert)) '. Median: ' num2str(median(expert))])
display(['Naive  horizon exp: ' num2str(mean(naive)) ' +- ' num2str(std(naive)) '. Median: ' num2str(median(naive))])
display(['Ranksum p=' num2str(p) ', z=', num2str(abs(stats.zval))])
display(['Cohen`s d=' num2str(d)])
[r, pCorr] = corr(asymptote(:), horizonExp, 'type', 'Spearman');
display(['Spearman corr=' num2str(r) ', p=' num2str(pCorr)])
text(30,20, ['R=' num2str(r,2) ' ***'])
nrep=1000;
rng(42);
for rep=1:nrep
    bootdiff(rep) = median(expert(randi(length(expert), length(expert), 1))) - median(naive(randi(length(naive), length(naive), 1))); 
end
display(['95% CI for the difference in medians: ' num2str(prctile(bootdiff, 2.5)) '--' num2str(prctile(bootdiff, 97.5))])

% [r,pCorr]=corr([asymptote(naiveGroup) - mean(asymptote(naiveGroup)) asymptote(expertGroup) - mean(asymptote(expertGroup))]', [horizonExp(naiveGroup) - mean(horizonExp(naiveGroup)); horizonExp(expertGroup) - mean(horizonExp(expertGroup))], 'type', 'Spearman')

% Slopes
subplot(236)
hold on
scatter(asymptote(expertGroup), slope(expertGroup), 25, expertColor, 'filled')
scatter(asymptote(naiveGroup), slope(naiveGroup), 25, naiveColor, 'filled')
ylim([0 7])
yticks(0:7)
xlim([20 100])
xticks(20:20:100)
xlabel('Time inside the path at the asymptote (%)')
ylabel('Initial slope (%/cm)')

% Compare naive and expert slopes
naive  = slope(naiveGroup);
expert = slope(expertGroup);
[p,~,stats] = ranksum(naive, expert);
d = abs(mean(naive)-mean(expert))/std([naive-mean(naive); expert-mean(expert)]);
disp(' ')
display(['Expert initial slope: ' num2str(mean(expert)) ' +- ' num2str(std(expert)) '. Median: ' num2str(median(expert))])
display(['Naive  initial slope: ' num2str(mean(naive)) ' +- ' num2str(std(naive)) '. Median: ' num2str(median(naive))])
display(['Ranksum p=' num2str(p) ', z=', num2str(abs(stats.zval))])
display(['Cohen`s d=' num2str(d)])
[r, pCorr] = corr(asymptote(:), slope, 'type', 'Spearman');
display(['Spearman corr=' num2str(r) ', p=' num2str(pCorr)])
text(30,5, ['R=' num2str(r,2) ' ***'])
nrep=1000;
rng(42);
for rep=1:nrep
    bootdiff(rep) = median(expert(randi(length(expert), length(expert), 1))) - median(naive(randi(length(naive), length(naive), 1))); 
end
display(['95% CI for the difference in medians: ' num2str(prctile(bootdiff, 2.5)) '--' num2str(prctile(bootdiff, 97.5))])

% Model prediction
E1=median(slope(naiveGroup)) * median(horizonChangepoint(naiveGroup)) + median(interceptChangepoint(naiveGroup));
E2=median(slope(expertGroup)) * median(horizonChangepoint(naiveGroup)) + median(interceptChangepoint(expertGroup));
E3=median(slope(expertGroup)) * median(horizonChangepoint(expertGroup)) + median(interceptChangepoint(expertGroup));
D = E3-E1;
F1 = (E2-E1)/(E3-E1);
F2 = (E3-E2)/(E3-E1);
display(['Median model estimates: ' num2str(E1) ', ' num2str(E2) ', ' num2str(D) ', ' num2str(F1) ', ' num2str(F2)])

E1=mean(slope(naiveGroup)) * mean(horizonChangepoint(naiveGroup)) + mean(interceptChangepoint(naiveGroup));
E2=mean(slope(expertGroup)) * mean(horizonChangepoint(naiveGroup)) + mean(interceptChangepoint(expertGroup));
E3=mean(slope(expertGroup)) * mean(horizonChangepoint(expertGroup)) + mean(interceptChangepoint(expertGroup));
D = E3-E1;
F1 = (E2-E1)/(E3-E1);
F2 = (E3-E2)/(E3-E1);
display(['Mean model estimates: ' num2str(E1) ', ' num2str(E2) ', ' num2str(D) ', ' num2str(F1) ', ' num2str(F2)])


% % Relationship with lags
% load 'lag.mat' % This file is made by figureTrajectories.m
% asymptoteLag = mean(lag(:,8:10), 2);
% 
% % Excluding 10 naive subjects with large asymptoteLag
% naive  = horizonChangepoint(naiveGroup & (asymptoteLag<10));
% expert = horizonChangepoint(expertGroup);
% [p,~,stats] = ranksum(naive, expert);
% d = abs(mean(naive)-mean(expert))/std([naive-mean(naive); expert-mean(expert)]);
% disp(' ')
% display(['Expert horizon: ' num2str(mean(expert)) ' +- ' num2str(std(expert)) '. Median: ' num2str(median(expert))])
% display(['Naive  horizon: ' num2str(mean(naive)) ' +- ' num2str(std(naive)) '. Median: ' num2str(median(naive))])
% display(['Ranksum p=' num2str(p) ', z=', num2str(abs(stats.zval))])
% display(['Cohen`s d=' num2str(d)])

