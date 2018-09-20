function linRegSpikes_opMD(xlFile)

[root, sep] = currComputer();

[~, sessionCellList, ~] = xlsread(xlFile);
sessionList = sessionCellList(2:end, 2);
cellList = sessionCellList(2:end, 1);

timeMax = 181000;
binSize = 30000;
timeBinEdges = [1000:binSize:timeMax];  %no trials shorter than 1s between outcome and CS on
tMax = length(timeBinEdges) - 1;

for i = 1:length(sessionList)

    %get spike information
    spikeStruct = [];
    spikeStruct = spikeProps_opMD(sessionList{i}, cellList{i});
    
    %load behavioral data
    [animalName, date] = strtok(sessionList{i}, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionList{i}(end), 'alpha')
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionList{i}(end) sep sessionList{i} '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionList{i} '_sessionData_behav.mat'];
    end

    if exist(sessionDataPath,'file')
        load(sessionDataPath)
    else
        [behSessionData, blockSwitch, blockSwitchL, blockSwitchR] = generateSessionData_operantMatchingDecoupled(sessionName);
    end

    %create arrays for choices and rewards
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    allChoice_R = double(allChoices == 1);
    allChoice_L = double(allChoices == -1);

    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = 1;

    %create binned outcome matrices
    rwdTmpMatx = NaN(tMax, length(responseInds));     %initialize matrices for number of response trials x number of time bins
    noRwdTmpMatx = NaN(tMax, length(responseInds));
    for j = 2:length(responseInds)          
        k = 1;
        %find time between "current" choice and previous rewards, up to timeMax in the past 
        timeTmpL = []; timeTmpR = []; nTimeTmpL = []; nTimeTmpR = [];
        while j-k > 0 & behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime < timeMax
            if behSessionData(responseInds(j-k)).rewardL == 1
                timeTmpL = [timeTmpL (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            if behSessionData(responseInds(j-k)).rewardR == 1
                timeTmpR = [timeTmpR (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            if behSessionData(responseInds(j-k)).rewardL == 0
                nTimeTmpL = [nTimeTmpL (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            if behSessionData(responseInds(j-k)).rewardR == 0
                nTtimeTmpR = [nTimeTmpR (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            k = k + 1;
        end
        %bin outcome times and use to fill matrices
        if ~isempty(timeTmpL)
            binnedRwds = discretize(timeTmpL,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedRwds == k)
                    rwdTmpMatx(k,j) = sum(binnedRwds == k);
                else
                    rwdTmpMatx(k,j) = 0;
                end
            end
        end
        if ~isempty(timeTmpR)
            binnedRwds = discretize(timeTmpR,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedRwds == k) & isnan(rwdTmpMatx(k,j))
                    rwdTmpMatx(k,j) = sum(binnedRwds == k);
                elseif ~isempty(binnedRwds == k) & ~isnan(rwdTmpMatx(k,j))
                    rwdTmpMatx(k,j) = rwdTmpMatx(k,j) + sum(binnedRwds == k);
                else
                    rwdTmpMatx(k,j) = 0;
                end
            end
        end
        if isempty(timeTmpL) & isempty(timeTmpR)
            rwdTmpMatx(:,j) = 0;
        end
        if ~isempty(nTimeTmpL)
            binnedNoRwds = discretize(nTimeTmpL,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedNoRwds == k)
                    noRwdTmpMatx(k,j) = sum(binnedNoRwds == k);
                else
                    noRwdTmpMatx(k,j) = 0;
                end
            end
        end
        if ~isempty(nTimeTmpR)
            binnedNoRwds = discretize(nTimeTmpR,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedNoRwds == k) & isnan(noRwdTmpMatx(k,j))
                    noRwdTmpMatx(k,j) = sum(binnedNoRwds == k);
                elseif ~isempty(binnedNoRwds == k) & ~isnan(noRwdTmpMatx(k,j))
                    noRwdTmpMatx(k,j) = noRwdTmpMatx(k,j) + sum(binnedNoRwds == k);
                else
                    noRwdTmpMatx(k,j) = 0;
                end
            end
        end
        if isempty(nTimeTmpL) & isempty(nTimeTmpR)
            noRwdTmpMatx(:,j) = 0;
        end
    end
    
    % linear regression models
    rwd_preCS = fitlm([rwdTmpMatx]', spikeStruct.preCScount);
    rwd_postCS = fitlm([rwdTmpMatx]', spikeStruct.postCScount);
    rwd_postCSrate = fitlm([rwdTmpMatx]', spikeStruct.maxCSrate);
    noRwd_preCS = fitlm([noRwdTmpMatx]', spikeStruct.preCScount);
    noRwd_postCS = fitlm([noRwdTmpMatx]', spikeStruct.postCScount);
    noRwd_postCSrate = fitlm([noRwdTmpMatx]', spikeStruct.maxCSrate);
    
    % plot beta coefficients from models
    figure; 
    subplot(1,3,1); title('pre CS'); hold on;
    relevInds = 2:tMax+1;
    coefVals = rwd_preCS.Coefficients.Estimate(relevInds);
    CIbands = coefCI(rwd_preCS);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    coefVals = noRwd_preCS.Coefficients.Estimate(relevInds);
    CIbands = coefCI(noRwd_preCS);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'-b', 'linewidth',2)
    xlabel('Outcome n seconds back')
    ylabel('\beta Coefficient')
    
    subplot(1,3,2); title('post CS'); hold on;
    relevInds = 2:tMax+1;
    coefVals = rwd_postCS.Coefficients.Estimate(relevInds);
    CIbands = coefCI(rwd_postCS);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    coefVals = noRwd_postCS.Coefficients.Estimate(relevInds);
    CIbands = coefCI(noRwd_postCS);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'-b', 'linewidth',2)
    xlabel('Outcome n seconds back')
    ylabel('\beta Coefficient')
    
    subplot(1,3,3); title('post CS rate'); hold on;
    relevInds = 2:tMax+1;
    coefVals = rwd_postCSrate.Coefficients.Estimate(relevInds);
    CIbands = coefCI(rwd_postCSrate);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    coefVals = noRwd_postCSrate.Coefficients.Estimate(relevInds);
    CIbands = coefCI(noRwd_postCSrate);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'-b', 'linewidth',2)
    xlabel('Outcome n seconds back')
    ylabel('\beta Coefficient')
    
   legend('Reward', 'No Reward')    
    xlim([0 tMax*binSize/1000 + binSize/1000])
    suptitle([sessionList{i} ' ' cellList{i}])
    set(gcf, 'Position', [-1596 335 1271 532])
    
end


