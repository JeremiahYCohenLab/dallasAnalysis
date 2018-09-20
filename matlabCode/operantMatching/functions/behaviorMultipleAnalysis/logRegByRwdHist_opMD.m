function [glm_rwdHistBin, lm_rwdHist, smoothWindows] = logRegByRwdHist_opMD(xlFile, animal, category)

[root, sep] = currComputer();

[weights, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

rwdHistBin = [];
rwdHist = [];
rwdsMatx = [];
choiceMatx = [];
combinedChangeChoice = [];
combinedAllChoices = [];
smoothWindows = [10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 180];
tMax = length(smoothWindows);

for i = 1: length(dayList)
    sessionName = dayList{i};
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end

    if exist(sessionDataPath,'file')
        load(sessionDataPath)
    else
        [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
    end
    
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    allChoice_R = double(allChoices == 1);
    allChoice_L = double(allChoices == -1);
    changeChoice = [false abs(diff(allChoices)) > 0];
    
    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allRewardsBin = zeros(1,length(allChoices));
    allRewardsBin(logical(allReward_R)) = 1;
    allRewardsBin(logical(allReward_L)) = 1;
    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = -1;
    rewardsList =  allRewards(find(allRewards~=0));
    
    %find all rwd times in terms of seconds
    choiceTimes = ceil(([behSessionData(responseInds).rewardTime] - behSessionData(1).CSon)/1000);     %baseline to start time and convert to s from ms
    rwdTimes = choiceTimes(logical(allRewardsBin));
    
    sessionTime = ceil((behSessionData(end).CSon + 3000 - behSessionData(1).CSon)/1000);     % find total session time and pad time for reward on last trial
    sessionRwdsBin = zeros(1,sessionTime);
    sessionRwdsBin(rwdTimes) = 1;
    sessionRwds = zeros(1,sessionTime);
    sessionRwds(rwdTimes(logical(rewardsList==1))) = 1;
    sessionRwds(rwdTimes(logical(rewardsList==-1))) = -1;

    sessionRwdsSmoothBin = [];
    sessionRwdsSmooth = [];
    for j = 1:length(smoothWindows)
        boxKern = ones(1,smoothWindows(j));                                       %smooth rewards over time
        sessionRwdsSmoothBinTmp = conv(sessionRwdsBin, boxKern);
        sessionRwdsSmoothTmp = conv(sessionRwds, boxKern);
        sessionRwdsSmoothBin(j,:) = sessionRwdsSmoothBinTmp(1:(end-(length(boxKern)-1)));
        sessionRwdsSmooth(j,:) = sessionRwdsSmoothTmp(1:(end-(length(boxKern)-1)));
    end
    
    choiceTimes = choiceTimes - 1;
    if choiceTimes(1) == 0
        rwdHistBinTmp = [zeros(length(smoothWindows),1) sessionRwdsSmoothBin(:,choiceTimes(2:end))];
        rwdHistTmp = [zeros(length(smoothWindows),1) sessionRwdsSmooth(:,choiceTimes(2:end))];
    else
        rwdHistBinTmp = sessionRwdsSmoothBin(:,choiceTimes);
        rwdHistTmp = sessionRwdsSmooth(:,choiceTimes);
    end
    
    %create matrices for choices and rwds tMax trials in the past
    rwdMatxTmp = [];
    choiceMatxTmp = [];
    for j = 1:tMax
        rwdMatxTmp(j,:) = [NaN(1,j) allRewards(1:end-j)];
        choiceMatxTmp(j,:) = [NaN(1,j) allChoices(1:end-j)];
    end
    
    rwdHistBin = [rwdHistBin NaN(length(smoothWindows),100) rwdHistBinTmp];
    rwdHist = [rwdHist NaN(length(smoothWindows),100) rwdHistTmp];
    rwdsMatx = [rwdsMatx NaN(tMax,100) rwdMatxTmp];
    choiceMatx = [choiceMatx NaN(tMax,100) choiceMatxTmp];
    combinedChangeChoice = [combinedChangeChoice NaN(1,100) changeChoice];
    combinedAllChoices = [combinedAllChoices NaN(1,100) allChoices];
    
end

%logistic regression models
glm_rwdHistBin = fitglm([rwdHistBin]', combinedChangeChoice,'distribution','binomial','link','logit'); 
rsq{1} = num2str(round(glm_rwdHistBin.Rsquared.Adjusted*100)/100);
lm_rwdHist = fitlm([rwdHist]', combinedAllChoices); 
rsq{2} = num2str(round(glm_rwdHistBin.Rsquared.Adjusted*100)/100);


%plot models
figure; hold on
relevInds = 2:length(smoothWindows)+1;
coefVals = glm_rwdHistBin.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rwdHistBin);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(smoothWindows+0.2,coefVals,errorL,errorU,'b','linewidth',2)

ylabel('\beta Coefficient')
xlim([0 (smoothWindows(end)+5)])
title([animal ' - Choice Switch on Reward History'])


%plot models
figure; hold on
relevInds = 2:length(smoothWindows)+1;
coefVals = lm_rwdHist.Coefficients.Estimate(relevInds);
CIbands = coefCI(lm_rwdHist);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(smoothWindows+0.2,coefVals,errorL,errorU,'b','linewidth',2)

ylabel('\beta Coefficient')
xlim([0 (smoothWindows(end)+5)])
title([animal ' - Choice on Side Reward History'])
