function [h] = rwdProbAtSwitch_opMD(xlFile, animal, category, timeWin, revForFlag, plotFlag)

if nargin < 6
    plotFlag = 0;
end
if nargin < 5
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

smoothKern = ones(1, timeWin);

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
    
    responseInds = find(~isnan([behSessionData.rewardTime]));
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL];
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    changeChoice = [false abs(diff(allChoices)) > 0];
    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allRewards = zeros(1,length(responseInds));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = 1;

    %% analyze reward rate as moving average
    choiceTimes = ceil(([behSessionData(responseInds).rewardTime] - behSessionData(1).CSon)/1000); %baseline to start time and convert to s from ms
    if choiceTimes(1) == 0
        choiceTimes(1) = 1;
    end
    rwdTimes_L = choiceTimes(logical(allReward_L));
    rwdTimes_R = choiceTimes(logical(allReward_R));

    sessionTime = ceil((behSessionData(end).CSon + 3000 - behSessionData(1).CSon)/1000);     % find total session time and pad time for reward on last trial
    sessionRwds = zeros(1,sessionTime);
    sessionRwds(rwdTimes_R) = 1;
    sessionRwds(rwdTimes_L) = -1;
    
    switchTimes = choiceTimes(changeChoice);
    
    smoothRwds = abs(conv(sessionRwds, smoothKern)) / timeWin;
    sessionAvg(i) =  mean(smoothRwds);
    rwdAtSwitch(i) = mean(smoothRwds(switchTimes));
    
    %% analyze in terms of logistic regression curve
    modelFolderLocation = [root animalName sep animalName 'sorted' sep];
    if exist(modelFolderLocation)
        modelFolder = dir(modelFolderLocation);
        modelInd = ~cellfun(@isempty,strfind({modelFolder.name},'glmTime.mat'));
        load([modelFolderLocation modelFolder(modelInd).name]);
    else
        if revForFlag == 1
            [~, ~, glm_all, t] = combineLogRegTime_opMD('probRevBeh.xlsx', animalName, 'all', revForFlag);
            mkdir(modelFolderLocation);
            save([modelFolderLocation animalName 'glmTime.mat'], 'glm_all', 't');
        else
            [~, ~, glm_all, t] = combineLogRegTime_opMD('goodBehDays.xlsx', animalName, 'preS', revForFlag);
            mkdir(modelFolderLocation);
            save([modelFolderLocation animalName 'glmTime.mat'], 'glm_all', 't');
        end
    end
    
    expFit = singleExpFit(glm_all.Coefficients.Estimate(2:t.tMax+1), [(1:t.tMax)*t.binSize/1000]);
    expConv = expFit.a*exp(-(1/expFit.b)*(1:t.timeMax/1000));
    expConv = expConv./sum(expConv);
    smoothRwdsLRM = conv(sessionRwds, expConv);
    smoothRwdsLRM = smoothRwdsLRM(1:(end-(length(expConv)-1)));
    sessionAvgLRM(i) = mean(abs(smoothRwdsLRM));
    rwdAtSwitchLRM(i) = mean(abs(smoothRwdsLRM(switchTimes)));
    
    
end

figure;
subplot(1,2,1); hold on; 
mdl = fitlm(sessionAvg, rwdAtSwitch);
plot(mdl)
legend({sprintf('Rsq = %i',mdl.Rsquared.Ordinary)})
plot([0 1], [0 1],'-k')
xlim([min(sessionAvg) max(sessionAvg)])
ylim([min(rwdAtSwitch) max(rwdAtSwitch)])
ylabel('rwds/s at switch'); xlabel('rwds/s session average')
title([sprintf('smooth window: %i, slope: %i', length(smoothKern), mdl.Coefficients.Estimate(2))])

subplot(1,2,2); hold on; 
mdl = fitlm(sessionAvgLRM, rwdAtSwitchLRM);
plot(mdl)
legend({sprintf('Rsq = %i',mdl.Rsquared.Ordinary)})
plot([0 1], [0 1],'-k')
xlim([min(sessionAvgLRM) max(sessionAvgLRM)])
ylim([min(rwdAtSwitchLRM) max(rwdAtSwitchLRM)])
ylabel('rwds/s at switch'); xlabel('rwds/s session average')
title([sprintf('slope: %i', mdl.Coefficients.Estimate(2))])
