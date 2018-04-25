clear all

rawDataFolder = '~/tmp/';
outputFileName = 'preprocessedData.mat';

% Expert subjects
E_subjects = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', ...
    's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'aa', 'ab', 'ac', 'ad', 'ae', 'af'};

% Naive subjects
N_subjects = {'an', 'bn', 'cn', 'dn', 'en', 'fn', 'gn', 'hn', 'in', 'jn', 'kn', 'ln', 'mn', 'nn', 'on', ...
    'pn', 'qn', 'rn', 'sn', 'tn', 'un', 'vn', 'wn', 'xn', 'yn', 'zn', 'aan', 'abn', 'acn', 'adn'};

expertGroup = logical([ones(length(E_subjects),1);  zeros(length(N_subjects), 1)]);
naiveGroup  = logical([zeros(length(E_subjects),1);  ones(length(N_subjects), 1)]);

% loop over subjects
for sub = 1:length(E_subjects)+length(N_subjects)
    fprintf('.')
    if sub <= length(E_subjects)
        d_sub = load([rawDataFolder 'Testing/Result_' E_subjects{sub} '_Testing.mat']);
    else
        d_sub = load([rawDataFolder 'Testing/Result_' N_subjects{sub-length(E_subjects)} '_Testing_NV.mat']);
    end
    
    % loop over trials
    for tr = 1:size(d_sub.Data,2)        
        x1 = d_sub.Data{tr}.expt_path(1:2:end,1);
        x2 = d_sub.Data{tr}.expt_path(1:2:end,3);
        c  = d_sub.Data{tr}.cursor_pos_path(1:2:end)';
        t  = d_sub.Data{tr}.timing(1:2:end)';
        
        midline = (x1+x2)/2;
        pathBorders = [x2 x1];
        pathBorders = bsxfun(@plus, midline, 50*50./bsxfun(@minus, pathBorders, midline));
        pathBorders = pathBorders * 1.1/40;
        
        expData{sub}{tr}.pathMidline = mean(pathBorders,2);
        expData{sub}{tr}.pathBorders = pathBorders;
        expData{sub}{tr}.cursorPosition = c * 1.1/40;
        expData{sub}{tr}.time = (t - t(1)) * 34.1;
        expData{sub}{tr}.searchlightLength = d_sub.Data{tr}.horizon;
    end
end
fprintf('\n')

% training data -- only expert subjects
for sub = 1:length(E_subjects)
    fprintf('.')
    
    %loop over training days
    for day = 1:5
        d_sub = load([rawDataFolder 'Training/Result_' E_subjects{sub} '_D' num2str(day) '_Training.mat']);
        
        % loop over trials
        for tr = 1:size(d_sub.Data,2)            
            x1 = d_sub.Data{tr}.expt_path(1:2:end,1);
            x2 = d_sub.Data{tr}.expt_path(1:2:end,3);
            c  = d_sub.Data{tr}.cursor_pos_path(1:2:end)';
            t  = d_sub.Data{tr}.timing(1:2:end)';
            
            midline = (x1+x2)/2;
            pathBorders = [x2 x1];
            pathBorders = bsxfun(@plus, midline, 50*50./bsxfun(@minus, pathBorders, midline));
            pathBorders = pathBorders * 1.1/40;
            
            trainingData{sub,day}{tr}.pathMidline = mean(pathBorders,2);
            trainingData{sub,day}{tr}.pathBorders = pathBorders;
            trainingData{sub,day}{tr}.cursorPosition = c * 1.1/40;
            trainingData{sub,day}{tr}.time = (t - t(1)) * 34.1;
            trainingData{sub,day}{tr}.searchlightLength = d_sub.Data{tr}.horizon;
        end
    end
end
fprintf('\n')

save(outputFileName, 'expData', 'expertGroup', 'naiveGroup', 'trainingData')
