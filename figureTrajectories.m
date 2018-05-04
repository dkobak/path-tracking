%% Analysis of time lags

clear all
load('preprocessedData.mat')

naiveColor = [228,26,28]/256;
expertColor = [55,126,184]/256;
cm = (10:10:100) / 100 * 30;

% sequence of searchlights (same for all subjects)
for trial = 1:length(expData{1})
    searchlightSeq(trial) = expData{1}{trial}.searchlightLength/10;
end

dt = 1/30 * 1000; % 30 Hz

for subject = 1:length(expData)
    for searchlight = 1:10       
        % combining all 3 trials with the same searchlight L together
        cursor = [];
        midline = [];        
        for trial = find(searchlightSeq == searchlight)
            cursor = [cursor; expData{subject}{trial}.cursorPosition];
            midline = [midline; expData{subject}{trial}.pathMidline];
        end
        
        % interpolating to increase resolution 10-fold
        cursor = interp(cursor, 10);
        midline = interp(midline, 10);
        
        % computing cross-correlation
        [r, lags] = xcorr(midline, cursor, 100, 'coeff');
        [~, indx] = max(r);
        lag(subject, searchlight) = -lags(indx) * dt / 10;
        
        % and now without combining
        trials = find(searchlightSeq == searchlight);
        for t = 1:3
            cursor =  expData{subject}{trials(t)}.cursorPosition;
            midline = expData{subject}{trials(t)}.pathMidline;
        
            cursor = interp(cursor, 10);
            midline = interp(midline, 10);
        
            [r, lags] = xcorr(midline, cursor, 100, 'coeff');
            [~, indx] = max(r);
            lag_trialwise(t, subject, searchlight) = -lags(indx) * dt / 10;
        end
    end
end

save('lag.mat', 'lag', 'lag_trialwise')

figure('Position', [100 100 1000 350]);
subplot(121)
hold on
plot(cm, lag(naiveGroup,:)', 'Color', naiveColor*.5 + [1 1 1]*.5)
plot(cm, lag(expertGroup,:)', 'Color', expertColor*.5 + [1 1 1]*.5)
plot(cm, mean(lag(naiveGroup,:)), '.-', 'Color', naiveColor, 'LineWidth', 3, 'MarkerSize', 30)
plot(cm, mean(lag(expertGroup,:)), '.-', 'Color', expertColor, 'LineWidth', 3, 'MarkerSize', 30)
axis([0 30 -100 300])
xticks([0 6 12 18 24 30])
xlabel('Searchlight length (cm)')
ylabel('Lag (ms)')

load accuracy.mat
asymptote = squeeze(mean(mean(accuracy(:,:,8:10),1),3));
asymptoteLag = mean(lag(:,8:10), 2);

subplot(122)
hold on
scatter(asymptote(expertGroup)', asymptoteLag(expertGroup)', 25, expertColor, 'filled')
scatter(asymptote(naiveGroup), asymptoteLag(naiveGroup), 25, naiveColor, 'filled')
axis([20 100 -100 300])
[r,p] = corr(asymptote', asymptoteLag, 'type', 'Spearman');
text(40,250,['R=' num2str(r,2) ' ***'])
disp(' ')
display(['Spearman corr=' num2str(r) ', p=' num2str(p)])

xlabel('Time inside the path at asymptote (%)')
ylabel('Asymptote lag (ms)') 


%% Analysis of path bends

figure('Position', [100 100 1600 600])

n = zeros(length(expData),10);
for searchlight = 1:10
    chunksPath = [];
    chunksPathU = [];
    chunksPathD = [];
    
    for trial = find(searchlightSeq == searchlight)
        borders = expData{1}{trial}.pathBorders;
        midline = expData{1}{trial}.pathMidline;

        window = 18;
        deep = 4.5;     
        % window=40, deep=7 yields around 10 segments per path
        % window=16, deep=4 yields around 15 per path
        % window=18, deep=4.5 yields around 13 per path
        
        minima{trial} = [];
        maxima{trial} = [];
       
        for pos = window : length(midline) - window
            before = midline((pos-window/2):(pos-1));
            after = midline((pos+1):(pos+window/2));
            if min([before; after]) >= midline(pos) ...
                    && max(before) > midline(pos) + deep ...
                    && max(after) > midline(pos) + deep
                minima{trial} = [minima{trial} pos];
            end
            
            if max([before; after]) <= midline(pos) ...
                    && min(before) < midline(pos) - deep ...
                    && min(after) < midline(pos) - deep
                maxima{trial} = [maxima{trial} pos];
            end
        end
        
        segmentsCount(trial) = length(minima{trial})+length(maxima{trial});
        
        for m = 1:length(minima{trial})
            win = minima{trial}(m)-window/2:minima{trial}(m)+window/2;
            chunksPath = [chunksPath; midline(win)' - midline(minima{trial}(m))];
            chunksPathU = [chunksPathU; borders(win,1)' - midline(minima{trial}(m))];
            chunksPathD = [chunksPathD; borders(win,2)' - midline(minima{trial}(m))];            
        end
        for m = 1:length(maxima{trial})
            win = maxima{trial}(m)-window/2:maxima{trial}(m)+window/2;
            chunksPath = [chunksPath; -midline(win)' + midline(maxima{trial}(m))];
            chunksPathU = [chunksPathU; -borders(win,2)' + midline(maxima{trial}(m))];
            chunksPathD = [chunksPathD; -borders(win,1)' + midline(maxima{trial}(m))];
        end
    end
    
    segmentsCountSearchl(searchlight) = size(chunksPath,1);
    
    subplot(2,5,searchlight)
    hold on
    x = expData{1}{1}.time(1:window+1);
    x = x - x(round(length(x)/2));
    plot(x, mean(chunksPathU), 'k', 'LineWidth', 2);
    plot(x, mean(chunksPathD), 'k', 'LineWidth', 2);
    
    axis([-10 10 -2 8])
    title(['L=' num2str(searchlight*10) '%'])
    interpFactor = 10;
    x = interp(x, interpFactor);
    x = x(1:end-interpFactor+1);
    
    for subject = 1:length(expData)
        chunksTraj = [];

        for trial = find(searchlight == searchlightSeq)
            midline = expData{1}{trial}.pathMidline;
            borders = expData{1}{trial}.pathBorders;
            cursor  = expData{subject}{trial}.cursorPosition;
            
            for m = 1:length(minima{trial})
                ind = minima{trial}(m) + (-window/2:window/2);
                chunksTraj = [chunksTraj; ...
                    cursor(ind)' - midline(minima{trial}(m))];

                acc = sum(cursor(ind)<borders(ind,1) & cursor(ind)>borders(ind,2)) ...
                    / length(ind) * 100;
                n(subject,searchlight) = n(subject,searchlight) + 1;
                accuracyChunks(subject,searchlight,n(subject,searchlight)) = acc;
            end
            for m = 1:length(maxima{trial})
                ind = maxima{trial}(m) + (-window/2:window/2);
                chunksTraj = [chunksTraj; ...
                    -cursor(ind)' + midline(maxima{trial}(m))];

                acc = sum(cursor(ind)<borders(ind,1) & cursor(ind)>borders(ind,2)) ...
                    / length(ind) * 100;
                n(subject,searchlight) = n(subject,searchlight) + 1;
                accuracyChunks(subject,searchlight,n(subject,searchlight)) = acc;
            end
        end
        
        bendPosStd(searchlight,subject) = mad(chunksTraj(:, window/2 + 1), 1);
        bendPosVer(searchlight,subject) = median(chunksTraj(:, window/2 + 1));
        
        chunksTrajInterp = [];
        for t = 1:size(chunksTraj,1)
            chint = interp(chunksTraj(t,:), interpFactor);
            chunksTrajInterp(t,:) = chint(1:end-interpFactor+1);
        end
        chunksTraj = chunksTrajInterp;
        
        [~, ind] = min(mean(chunksTraj));
        bendPoints(subject,:) = [x(ind) mean(chunksTraj(:,ind))];
        bendPosStdAlt(searchlight,subject) = mad(chunksTraj(:, ind), 1);
        bendPosVerAlt(searchlight,subject) = mean(chunksTraj(:,ind));
        
        if naiveGroup(subject)
            plot(x, mean(chunksTraj), 'Color', [1 .5 .5])%, 'LineWidth', 2)
            plot(x(ind), mean(chunksTraj(:,ind)), 'r.', 'MarkerSize', 20)
        else
            plot(x, mean(chunksTraj), 'Color', [.4 .7 .4])%, 'LineWidth', 2)
            plot(x(ind), mean(chunksTraj(:,ind)), '.', 'Color', [0 .6 0], 'MarkerSize', 20)           
        end
    end
    
    % code for drawing coverage contours
    
    [bandwidth,density,X,Y]=kde2d(bendPoints(naiveGroup,:));
    %figure; surf(X,Y,density,'LineStyle','none')
    v = .01:.001:.2;
    for i=1:length(v)    
        integral(i) = (X(1,2)-X(1,1))*(Y(2,1)-Y(1,1))*sum(density(density>v(i)));
    end
    vCutoff = v(find(integral<.75, 1));  
    contour(X,Y,density, [vCutoff vCutoff], 'Color', [1 0 0], 'LineWidth', 2)
    [bandwidth,density,X,Y]=kde2d(bendPoints(expertGroup,:));
    %figure; surf(X,Y,density,'LineStyle','none')
    v = .01:.001:.2;
    for i=1:length(v)    
        integral(i) = (X(1,2)-X(1,1))*(Y(2,1)-Y(1,1))*sum(density(density>v(i)));
    end
    vCutoff = v(find(integral<.75, 1));    
    contour(X,Y,density, [vCutoff vCutoff], 'Color', [0 .6 0], 'LineWidth', 2)
    
    pause(0.1)
end

display(' ')
display(['Segments per path: ' num2str(mean(segmentsCount)) ' +- ' num2str(std(segmentsCount))])
display(['Segments per path: ' num2str(mean(segmentsCountSearchl)) ' +- ' num2str(std(segmentsCountSearchl))])


