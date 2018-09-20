function [glm_rwdLick] = combineLogRegLickLat_opM(xlFile, animal, category)

%make flexible for M and MD, run on lesioned animals

[root, sep] = currComputer();

[weights, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end
combinedRewardsMatx = [];
combinedLickLat = [];
tMax = 10;


for i = 1: length(dayList)
    sessionName = dayList{i};
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName(end) sep sessionName '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sessionName '_sessionData_behav.mat'];
    end

    if exist(sessionDataPath,'file')
        load(sessionDataPath)
    else
        [sessionData, ~] = generateSessionData_operantMatching(sessionName);
    end
    
    %%generate reward matrix for tMax trials
    responseInds = find(~isnan([sessionData.rewardTime])); % find CS+ trials with a response in the lick window
    allReward_R = [sessionData(responseInds).rewardR]; 
    allReward_L = [sessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(sessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = 1;

    
    rwdMatxTmp = [];
    for j = 1:tMax
        rwdMatxTmp(j,:) = [NaN(1,j) allRewards(1:end-j)];
    end

    combinedRewardsMatx = [combinedRewardsMatx NaN(tMax,100) rwdMatxTmp(:,1:end-1)];
    
    
    %% determine and plot lick latency distributions for each spout
    lickLat = [sessionData(responseInds).rewardTime] - [sessionData(responseInds).CSon];
    indsR = find(allChoices == 1);
    indsL = find(allChoices == -1);
    lickLat_R = zscore(lickLat(indsR));
    lickLat_L = zscore(lickLat(indsL));
    lickLat = NaN(1, length(allChoices));
    lickLat(indsR) = lickLat_R;
    lickLat(indsL) = lickLat_L;
    
    combinedLickLat = [combinedLickLat NaN(1,100) lickLat(2:end)];

end

%linear regression model
glm_rwdLick = fitlm([combinedRewardsMatx]', combinedLickLat);

figure; hold on;
relevInds = 2:tMax+1;
coefVals = glm_rwdLick.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rwdLick);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color', [1 0 1],'linewidth',2)

xlabel('Reward n Trials Back')
ylabel('\beta Coefficient')
xlim([0.5 tMax+0.5])

suptitle([animal ' ' category])