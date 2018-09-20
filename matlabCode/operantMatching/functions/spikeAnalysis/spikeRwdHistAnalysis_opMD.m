function [tbl] = spikeRwdHistAnalysis_opMD(cellFile, category, plotFlag)

if nargin < 3
    plotFlag = 0;
end

[root, sep] = currComputer();

[revForFlag, cellList, ~] = xlsread(cellFile, category);
cellList = cellList(2:end,:);

maxXCF = NaN(length(cellList), 1);
maxLag = NaN(length(cellList), 1);
rowNames = cell(length(cellList), 1);

binEdges = [0:0.2:1];


for i = 1:length(cellList)
    sessionName = cellList{i,2};
    [animalName] = strtok(sessionName, 'd');
    animalName = animalName(2:end);

    if isstrprop(sessionName(end), 'alpha')
        sortedFolderLocation = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
    else
        sortedFolderLocation = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
    end
    sortedFolder = dir(sortedFolderLocation);


%     if revForflag{i} == 1
%         if any(~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'))) == 1
%             sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'));
%             load([sortedFolderLocation sortedFolder(sessionDataInd).name])
%         else
%             [sessionData] = generateSessionData_intan_operantMatching(sessionName);
%         end
%     else
        if any(~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'))) == 1
            sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'));
            load([sortedFolderLocation sortedFolder(sessionDataInd).name])
        else
            [sessionData] = generateSessionData_nL_operantMatching(sessionName);
        end
%    end
    [s] = behAnalysisNoPlot_opMD(sessionName, 'revForFlag', revForFlag(i));

%     if ~isempty(p.Results.trialList)
%         s.responseInds = s.responseInds(p.Results.trialList);
%         s.sessionRwds = s.sessionRwds(p.Results.trialList);
%         rwdHx = rwdHx(p.Results.trialList);
%     end
    
    %% smooth rewards over time
    modelFolderLocation = [root animalName sep animalName 'sorted' sep];
    if exist(modelFolderLocation)
        modelFolder = dir(modelFolderLocation);
        modelInd = ~cellfun(@isempty,strfind({modelFolder.name},'glmTime.mat'));
        load([modelFolderLocation modelFolder(modelInd).name]);
    else
        if revForFlag(i) == 1
            [~, ~, glm_all, t] = combineLogRegTime_opMD('probRevBeh.xlsx', animalName, 'all', revForFlag(i));
            mkdir(modelFolderLocation);
            save([modelFolderLocation animalName 'glmTime.mat'], 'glm_all', 't');
        else
            [~, ~, glm_all, t] = combineLogRegTime_opMD('goodBehDays.xlsx', animalName, 'all', revForFlag(i));
            mkdir(modelFolderLocation);
            save([modelFolderLocation animalName 'glmTime.mat'], 'glm_all', 't');
        end
    end
    
    expFit = singleExpFit(glm_all.Coefficients.Estimate(2:t.tMax+1), [(1:t.tMax)*t.binSize/1000]);
    expConv = expFit.a*exp(-(1/expFit.b)*(1:t.timeMax/1000));
    expConv = expConv./sum(expConv);
    sessionRwdsSmooth = conv(s.sessionRwds, expConv);
    sessionRwdsSmooth = sessionRwdsSmooth(1:(end-(length(expConv)-1)));

    
    %% smooth spikes over the session
    fieldNames = fieldnames(sessionData);
    cellInd = ~cellfun(@isempty,strfind(fieldNames, cellList{i,1}));
    cellInd = find(cellInd == 1);
    cellInd = cellInd - 9;
    sessionSpikes = sessionData(cellInd).allSpikes - sessionData(1).CSon;
    sessionSpikeInd = [sessionSpikes(sessionSpikes > 0)] / 1000;
    sessionSpikeCounts = histcounts(sessionSpikeInd, [0:length(s.sessionRwds)]);    
    
    [xcf, lags] = crosscorr(sessionSpikeCounts, sessionRwdsSmooth);
    if max(abs(xcf)) > max(xcf)
        [maxXCF(i), maxInd] = min(xcf);
    else
        [maxXCF(i), maxInd] = max(xcf);
    end
    maxLag(i) = lags(maxInd);
    rowNames{i} = [cellList{i,2} '_' cellList{i,1}];
    
    if maxXCF(i) > 0.2 | maxXCF(i) < -0.2
        sessionRwdsSmoothNorm = sessionRwdsSmooth/max(sessionRwdsSmooth);
        rwdHxSpikes =[]; rwdHxSpikesSD = [];
        for j = 1:length(binEdges)-1
            tmp = sessionSpikeCounts(sessionRwdsSmoothNorm >= binEdges(j) & sessionRwdsSmoothNorm < binEdges(j+1));
            rwdHxSpikes(j) = mean(tmp);
            rwdHxSpikesSD(j) = std(tmp) / sqrt(length(tmp));
        end
        meow = figure;
        subplot(1,2,1)
        scatter(sessionSpikeCounts, sessionRwdsSmooth);
        xlabel('spike counts'); ylabel('smoothed reward history')
        subplot(1,2,2)
        errorbar([binEdges(1:end-1)+diff(binEdges)/2], rwdHxSpikes, rwdHxSpikesSD)
        xlabel('smoothed and normalized reward history'); ylabel('spike counts')
        suptitle([cellList{i,2} cellList{i,1}])
        saveFigurePDF(meow,['C:\Users\cooper_PC\Desktop\analysis\rwdHxSpikeCounts\rwdHxSpikeCounts' num2str(i)]);
    end
        
end
 
tbl = table(maxXCF, maxLag, 'VariableNames', {'maxXCF' 'maxLag'}, 'RowNames', rowNames); 

if plotFlag
    figure; 
    subplot(1,2,1)
    histogram(tbl.maxXCF, 'Normalization', 'Probability', 'FaceColor', 'c')
    ylabel('probability')
    xlabel('max cross correlation coeff')

    subplot(1,2,2)
    histogram(tbl.maxLag, 'Normalization', 'Probability', 'FaceColor', 'm')
    ylabel('probability')
    xlabel('max cross correlation lag')
end
