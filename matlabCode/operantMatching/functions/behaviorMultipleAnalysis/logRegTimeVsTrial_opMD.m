function [glm_rwd, glm_noRwd, glm_all, binSize, timeMax] = logRegTimeVsTrial_opMD(xlFile, animal, category)

%determine root for file location
[root, sep] = currComputer();

%get values for exponential decay filters
[~, ~, glm_rwd] = combineLogReg_opMD(xlFile, animal, category); close;      %for trial-wise
expFit = singleExpFit(glm_rwd.Coefficients.Estimate(2:end));
expConvTrial = expFit.a*exp(-(1/expFit.b)*(1:10));
expConvTrial = expConvTrial./sum(expConvTrial);

[glm_rwdTime, ~, ~, binSize, timeMax] = logRegTime_opMD(xlFile, animal, category); close;       %for time-wise
xVals = [1:binSize/1000:timeMax/1000];
xVals = xVals(1:end-1) + diff(xVals)/2;
expFitTime = singleExpFit(glm_rwdTime.Coefficients.Estimate(2:end), xVals');
expConvTime = expFitTime.a*exp(-(1/expFitTime.b)*(1:timeMax/1000));



%import behavior session titles for desired category
[weights, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

%initialize matrices
rwdHx = []; 
rwdHxTimeChoice = [];
combinedChangeChoice = [];

%loop for each session in the list
for i = 1: length(dayList)              
    sessionName = dayList{i};                       %extract relevant info from session title
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')         %define appropriate data path
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end

    if exist(sessionDataPath,'file')        %load preprocessed struct if there is one
        load(sessionDataPath)
    else                                    %otherwise generate the struct
        [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
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
    changeChoice = [false abs(diff(allChoices)) > 0];
    
    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = 1;
    allRewardsBinary = allRewards;                      %make all rewards have the same value
    allRewardsBinary(find(allRewards==-1)) = 1;
    rewardsList =  allRewards(find(allRewards~=0));
    
    %convolve with trial-wise filter
    rwdHxTmp = conv(allRewardsBinary,expConvTrial);              %convolve with exponential decay to give weighted moving average
    rwdHxTmp = rwdHxTmp(1:end-(length(expConvTrial)-1));  
    
   
    %find all rwd times in terms of seconds
    choiceTimes = ceil(([behSessionData(responseInds).rewardTime] - behSessionData(1).CSon)/1000);     %baseline to start time and convert to s from ms
    rwdTimes = choiceTimes(logical(allRewardsBinary));

    sessionTime = ceil((behSessionData(end).CSon + 3000 - behSessionData(1).CSon)/1000);     % find total session time and pad time for reward on last trial
    sessionRwdsBin = zeros(1,sessionTime);
    sessionRwdsBin(rwdTimes) = 1;
    sessionRwds = zeros(1,sessionTime);
    sessionRwds(rwdTimes(logical(rewardsList==1))) = 1;
    sessionRwds(rwdTimes(logical(rewardsList==-1))) = 1;
    
    %convolve with time-wise filter
    rwdHxTimeTmp = conv(sessionRwds,expConvTime);              %convolve with exponential decay to give weighted moving average
    rwdHxTimeTmp = rwdHxTimeTmp(1:end-(length(expConvTime)-1)); 
    
    %take rwdHx at time of choice
    choiceTimes = choiceTimes - 1;
    if choiceTimes(1) == 0
        rwdHxTimeChoiceTmp = [0 rwdHxTimeTmp(choiceTimes(2:end))];
    else
        rwdHxTimeChoiceTmp = rwdHxTimeTmp(:,choiceTimes);
    end
    
    
    %combine with other sessions
    
    rwdHx = [rwdHx NaN(1,100) rwdHxTmp];
    rwdHxTimeChoice = [rwdHxTimeChoice NaN(1,100) rwdHxTimeChoiceTmp];
    combinedChangeChoice = [combinedChangeChoice NaN(1,100) changeChoice];
    
end

%logistic regression models
glm_all = fitglm([rwdHx' rwdHxTimeChoice'], combinedChangeChoice,'distribution','binomial','link','logit'); rsq{2} = num2str(round(glm_rwd.Rsquared.Adjusted*100)/100);

figure; hold on;
coefVals = glm_rwd.Coefficients.Estimate(2);
CIbands = coefCI(glm_all);
errorL = abs(coefVals - CIbands(2,1));
errorU = abs(coefVals - CIbands(2,2));
errorbar([1],coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)

coefVals = glm_rwd.Coefficients.Estimate(3);
CIbands = coefCI(glm_all);
errorL = abs(coefVals - CIbands(3,1));
errorU = abs(coefVals - CIbands(3,2));
errorbar([2],coefVals,errorL,errorU,'b','linewidth',2)

ylabel('\beta Coefficient')
xlim([0 3])
xticks([1 2])
xticklabels({'trial', 'time'})
title([animal ' ' category])