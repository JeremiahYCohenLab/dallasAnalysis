function [glm_rwdNoRwd, tMax] = combineLogReg_opMD(xlFile, animal, category, revForFlag, plotFlag)

if nargin < 5
    plotFlag = 0;
end
if nargin < 4
    revForFlag = 0;
end

[root, sep] = currComputer();

[weights, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end
rwdRateMatx = [];
combinedChoicesMatx = []; 
combinedRewardsMatx = [];
combinedNoRewardsMatx = [];
combinedTimesMatx = [];
combinedAllChoice_R = [];
tMax = 12;
 

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
        if revForFlag
            behSessionData = sessionData;
        end
    elseif revForFlag                                    %otherwise generate the struct
        [behSessionData, ~] = generateSessionData_operantMatching(sessionName);
    else
        [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
    end
    
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    omitInds = isnan([behSessionData.rewardTime]); 
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    
    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allChoice_R = double(allChoices == 1);
    allChoice_L = double(allChoices == -1);

    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = -1;
    
    allNoRewards = allChoices;
    allNoRewards(logical(allReward_R)) = 0;
    allNoRewards(logical(allReward_L)) = 0;
    
    outcomeTimes = [behSessionData(responseInds).rewardTime] - behSessionData(responseInds(1)).rewardTime;
    outcomeTimes = [diff(outcomeTimes) NaN];
    
    rwdMatxTmp = [];
    choiceMatxTmp = [];
    noRwdMatxTmp = [];
    for j = 1:tMax
        rwdMatxTmp(j,:) = [NaN(1,j) allRewards(1:end-j)];
        choiceMatxTmp(j,:) = [NaN(1,j) allChoices(1:end-j)];
        noRwdMatxTmp(j,:) = [NaN(1,j) allNoRewards(1:end-j)];
    end

    timeTmp = NaN(tMax,length(allRewards)); 
    for j = 1:tMax
        for k = 1:length(outcomeTimes)-j
            timeTmp(j,k+j) = sum(outcomeTimes(k:k+j-1));
        end
    end
    
%    allRewards(allRewards == -1) = 1;
    rwdsTmp = NaN(tMax,length(allRewards)); 
    for j = 1:tMax
        for k = 1:length(outcomeTimes)-j
            rwdsTmp(j,k+j) = sum(allRewards(k:k+j-1));
        end
    end
    
    rwdRateMatx = [rwdRateMatx NaN(tMax, 100) (rwdsTmp ./ timeTmp)];
    combinedRewardsMatx = [combinedRewardsMatx NaN(tMax,100) rwdMatxTmp];
    combinedNoRewardsMatx = [combinedNoRewardsMatx NaN(tMax,100) noRwdMatxTmp];
    combinedChoicesMatx = [combinedChoicesMatx NaN(tMax,100) choiceMatxTmp];
    combinedTimesMatx = [combinedTimesMatx NaN(tMax, 100) timeTmp];
    combinedAllChoice_R = [combinedAllChoice_R NaN(1,100) allChoice_R];
end


%logistic regression models
glm_rwd = fitglm([combinedRewardsMatx]', combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{1} = num2str(round(glm_rwd.Rsquared.Adjusted*100)/100);
glm_choice = fitglm([combinedChoicesMatx]', combinedAllChoice_R, 'distribution','binomial','link','logit'); rsq{3} = num2str(round(glm_choice.Rsquared.Adjusted*100)/100);
glm_rwdANDchoice = fitglm([combinedRewardsMatx; combinedChoicesMatx]', combinedAllChoice_R, 'distribution','binomial','link','logit'); rsq{2} = num2str(round(glm_rwdANDchoice.Rsquared.Adjusted*100)/100);
glm_time = fitglm([combinedTimesMatx]', combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{4} = num2str(round(glm_time.Rsquared.Adjusted*100)/100);
glm_rwdANDtime = fitglm([combinedRewardsMatx; combinedTimesMatx]', combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{5} = num2str(round(glm_rwdANDtime.Rsquared.Adjusted*100)/100);
glm_rwdRate = fitglm([rwdRateMatx]', combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{6} = num2str(round(glm_rwd.Rsquared.Adjusted*100)/100);
glm_rwdNoRwd = fitglm([combinedRewardsMatx; combinedNoRewardsMatx]', combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{7} = num2str(round(glm_rwdNoRwd.Rsquared.Adjusted*100)/100);
% glm_noRwdANDchoice = fitglm([combinedNoRewardsMatx; combinedChoicesMatx]', combinedAllChoice_R, 'distribution','binomial','link','logit'); rsq{8} = num2str(round(glm_noRwdANDchoice.Rsquared.Adjusted*100)/100);


if plotFlag
    figure; hold on;
    relevInds = 2:tMax+1;
    coefVals = glm_rwdNoRwd.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdNoRwd);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)

    relevInds = tMax+2:length(glm_rwdNoRwd.Coefficients.Estimate);
    coefVals = glm_rwdNoRwd.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdNoRwd);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'b','linewidth',2)

    xlabel('Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])
    title([animal ' ' category])
end
