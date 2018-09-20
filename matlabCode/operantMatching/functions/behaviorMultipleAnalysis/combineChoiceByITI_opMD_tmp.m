 function [x, probSwitchNoRwd, probStayRwd] = combineChoiceByITI_opMD_tmp(xlFile, animal, category, revForFlag)

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

ITIthresh = [0 10 100];

ITIcomb = [];
changeChoiceComb = [];
allRewardsComb = [];

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

    if revForFlag
        if exist(sessionDataPath,'file')
            load(sessionDataPath)
            behSessionData = sessionData;
        else
            [behSessionData, blockSwitch, blockProbs] = generateSessionData_behav_operantMatching(sessionName);
        end
    else
        if exist(sessionDataPath,'file')
            load(sessionDataPath)
        else
            [behSessionData, blockSwitch, blockSwitchL, blockSwitchR] = generateSessionData_operantMatchingDecoupled(sessionName);
        end
    end
    
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    
    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;

    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = 1;
    allRewards = allRewards(1:end-1);
    
    ITI = nan(1,length(responseInds) - 1); 
    for j = 2:length(responseInds)
        ITI(j-1) = behSessionData(responseInds(j)).CSon - behSessionData(responseInds(j-1)).rewardTime;
    end
    ITIcomb = [ITIcomb ITI/1000];
        
    %combine arrays
    changeChoice = [abs(diff(allChoices)) > 0];
    changeChoiceComb = [changeChoiceComb changeChoice];
    allRewardsComb = [allRewardsComb allRewards];  
    
end

for i = 1:length(ITIthresh)-1
    indsTmp = find(ITIcomb > ITIthresh(i) & ITIcomb < ITIthresh(i+1));
    changeChoiceTmp = changeChoiceComb(indsTmp);
    allRewardsTmp = allRewardsComb(indsTmp);
    probSwitchNoRwd(i) = sum(changeChoiceTmp(allRewardsTmp==0))/sum(allRewardsTmp==0);
    probStayRwd(i) = 1 - (sum(changeChoiceTmp(allRewardsTmp==1))/sum(allRewardsTmp==1));
end

figure; hold on;
x = ITIthresh(1:end-1) + diff(ITIthresh)/2;
plot(x, probSwitchNoRwd, '-m', 'linewidth', 2)
plot(x, probStayRwd, '-c', 'linewidth', 2)
legend('P(stay|rwd)', 'P(switch|no rwd)')
title([animal ' ' category])




